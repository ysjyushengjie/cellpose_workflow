version 1.0

# Cellpose 工作流：用于细胞分割和单细胞矩阵生成以及报告生成
workflow single_cell_analysis {
    input {
        # 输入的 TIF 图像文件 (HE染色图像)
        File input_tif
        # 输入的 GEM 文件
        File input_gem
        # 芯片 ID
        String chip_id
        # Cellpose 模型路径
        String model_path
        # Docker 镜像地址
        String docker_image
    }

    # GEM转换为TIFF
    call gem2tiff {
        input:
            input_gem = input_gem,
            chip_id = chip_id,
            docker_image = docker_image
    }

    # 细胞分割
    call segment_cell {
        input:
            input_tif = input_tif,
            model_path = model_path,
            chip_id = chip_id,
            docker_image = docker_image
    }

    # 生成单细胞矩阵
    call generate_matrix {
        input:
            input_gem = input_gem,
            chip_id = chip_id,
            mask_path = segment_cell.mask_output,
            docker_image = docker_image
    }

    # 生成HTML报告
    call generate_report {
        input:
            chip_id = chip_id,
            h5ad_matrix = generate_matrix.h5ad_output,
            gem_tiff = gem2tiff.tiff_output,
            docker_image = docker_image
    }

    output {
        # GEM可视化文件
        File gem_tiff = gem2tiff.tiff_output
        File gem_png = gem2tiff.png_output
        # 细胞掩码文件
        File mask_file = segment_cell.mask_output
        # 单细胞矩阵
        File h5ad_matrix = generate_matrix.h5ad_output
        File cell_bin_gef = generate_matrix.cell_bin_gef
        # HTML报告
        File html_report = generate_report.report_html
        File experiment_data_json = generate_report.experiment_data_json
    }
}

# GEM转TIFF任务
task gem2tiff {
    input {
        File input_gem
        String chip_id
        String docker_image
    }

    command <<<
        cat > 01_gem2tiff.py << 'EOL'
import numpy as np
from PIL import Image
import spateo as st
import matplotlib.pyplot as plt
import sys

def convert_gem2tiff(input_file, chip_id):
    adata = st.io.read_bgi_agg(input_file)
    img = adata.X.toarray()
    img[img > 10] = 10
    img = img / 10 * 255 # Convert to grayscale image
    
    plot_data = Image.fromarray(adata.X.toarray().T)

    if img.shape == adata.shape:
        img = Image.fromarray(img.T)
        img = img.convert("L")
        img.save(f'./{chip_id}_bin1_beforePS.tiff')
    else:
        raise ValueError(f"Shape mismatch: img.shape is {img.shape}, but adata.shape is {adata.shape}.\nPlease check your input data")
        
    return plot_data
        
def plot_chip(plot_data, chip_id):
    plt.figure(figsize=(6, 8))
    plt.imshow(plot_data, cmap='viridis', vmax='1')
    plt.axis('off')
    plt.savefig(f'./{chip_id}_bin1_beforePS.png', format='png')
    
    
if __name__ == "__main__":
    # parameter infomation
    input_file = sys.argv[1] # for example: '/data/input/Files/ResultData/Workflow/T202412180002227/04.tissuecut/Y00724H9.gem.gz'
    chip_id = sys.argv[2] # Y00724H9
    plot_data = convert_gem2tiff(input_file, chip_id)
    plot_chip(plot_data, chip_id)
EOL

        python 01_gem2tiff.py ~{input_gem} ~{chip_id}
    >>>

    output {
        File tiff_output = "${chip_id}_bin1_beforePS.tiff"
        File png_output = "${chip_id}_bin1_beforePS.png"
    }

    runtime {
        docker: docker_image
        memory: "16 GB"
        cpu: 4
    }
}

# 细胞分割任务
task segment_cell {
    input {
        File input_tif
        String model_path
        String chip_id
        String docker_image
    }

    command <<<
        cat > 02_based_cellpose_segment_cell.py << 'EOL'
from cellpose import models, io
import matplotlib.pyplot as plt
import sys


def cell_segment(registered_file, chip_id, model_path='/data/work/ascites/model/models/cyto3'):
    # read registered file
    img = io.imread(registered_file)
    # load pre-trained model
    model = models.CellposeModel(gpu=True, pretrained_model=model_path)
    # run pre-trained model
    masks, flows, styles = model.eval(img, diameter=None, channels=[0,3])

    io.masks_flows_to_seg(img, masks, flows, registered_file, styles, channels=[0,3])
    io.save_masks(img, masks, flows, f'{chip_id}_cellpose', png=True, tif=True, channels=[0,3])


if __name__ == "__main__":
    # regesited tiff, for example '/data/work/ascites/aligned_picture/Y00722A9_HE_aligned_OYK.tif'
    input_file = sys.argv[1]
    model_path = sys.argv[2]
    chip_id = sys.argv[3]
    
    cell_segment(input_file, chip_id, model_path)
EOL

        python 02_based_cellpose_segment_cell.py ~{input_tif} ~{model_path} ~{chip_id}
    >>>

    output {
        File mask_output = "${chip_id}_cellpose_cp_masks.png"
    }

    runtime {
        docker: docker_image
        memory: "16 GB"
        cpu: 4
        gpu: true
    }
}

# 生成单细胞矩阵任务
task generate_matrix {
    input {
        File input_gem
        String chip_id
        File mask_path
        String docker_image
    }

    command <<<
        cat > 03_cell_correct_and_generate_matrix.py << 'EOL'
import stereo as st
from stereo.tools.cell_correct import cell_correct
from stereo.tools.cell_cut import CellCut
import gzip
import sys


def gem_file_check(gem_path, chip_id):
    with gzip.open(gem_path, 'rt', encoding='utf-8') as file:
        lines = [next(file) for _ in range(8)]
        rest_of_info = file.read() 
        
        print(lines)
        
    modified_status = False 
    if lines[-1].startswith("#OffsetY=") and lines[-1] != "#OffsetY=0\n":
        lines[-1] = "#OffsetY=0\n"
        modified_status = True
    if lines[-2].startswith("#OffsetX=") and lines[-2] != "#OffsetX=0\n":
        lines[-2] = "#OffsetX=0\n"
        modified_status = True
    
    if modified_status:
        print(f'gem.gz file has been modified: {modified_status}')

        with gzip.open(f'./{chip_id}_modified.gem.gz', 'wt', encoding='utf-8') as file:
            file.writelines(lines)
            file.writelines(rest_of_info)
    
    return modified_status

                
def cellcut(gem_path, mask_path, chip_id, modified_status):
    cgef_out_dir = "./"
    cc = CellCut(cgef_out_dir=cgef_out_dir)
    if modified_status:
        gem_path = f'./{chip_id}_modified.gem.gz'
        out_path = cc.cell_cut(gem_path=gem_path, mask_path=mask_path)
    else:
        out_path = cc.cell_cut(gem_path=gem_path, mask_path=mask_path)

def cellcorrect(gem_path, mask_path, chip_id, modified_status):
    out_dir = "./"
    if modified_status:
        gem_path = f'./{chip_id}_modified.gem.gz'
        data = cell_correct(out_dir=out_dir, gem_path=gem_path, mask_path=mask_path, only_save_result=False, fast='FAST')
    else:
        data = cell_correct(out_dir=out_dir, gem_path=gem_path, mask_path=mask_path, only_save_result=False, fast='FAST')
    

def gef2matrix(chip_id):
    data = st.io.read_gef(f'./{chip_id}.adjusted.cellbin.gef', bin_type='cell_bins')
    st.io.stereo_to_anndata(data, flavor='scanpy',output=f'./{chip_id}_cellpose.h5ad')
    
    
if __name__ == "__main__":
    # like Y00722CC.gem.gz
    gem_path = sys.argv[1]
    # like Y00722CC
    chip_id = sys.argv[2]
    # the mask file is after cell segmentation
    mask_path = sys.argv[3]
    
    print('gem file check---------started')
    modified_status = gem_file_check(gem_path, chip_id)
    print('gem file check---------finished')
    print('cellcut---------started')
    cellcut(gem_path, mask_path, chip_id, modified_status)
    print('cellcut---------finished')
    print('cellcorrect----------started')
    cellcorrect(gem_path, mask_path, chip_id, modified_status)
    print('cellcorrect----------finished')
    print('generate singlecell matrix-------------started')
    gef2matrix(chip_id)
    print('generate singlecell matrix-------------finished')
EOL

        python 03_cell_correct_and_generate_matrix.py ~{input_gem} ~{chip_id} ~{mask_path}
    >>>

    output {
        File h5ad_output = "${chip_id}_cellpose.h5ad"
        File cell_bin_gef = "${chip_id}.adjusted.cellbin.gef"
    }

    runtime {
        docker: docker_image
        memory: "16 GB"
        cpu: 4
    }
}

# 生成HTML报告任务
task generate_report {
    input {
        String chip_id
        File h5ad_matrix
        File gem_tiff
        String docker_image
    }

    command <<<
        # 创建目录结构
        mkdir -p img interactive
        cp ~{gem_tiff} img/gem_visualization.png

        # 创建报告模板文件
        cat > report_template.md << 'EOL'
# 病理诊断模型报告

## 病人信息

| 病人姓名   | {{patient_name}} | 实验类型   | 芯片单细胞 |
|------------|-------------------|------------|------------|
| 病理号     | {{patient_id}}    | 样本类型   | 腹腔积液   |
| 送检日期   | {{sample_date}}   | 报告日期   | {{report_date}} |
| 临床病理诊断 | {{clinical_diagnosis}} | | |

## 实验信息

![实验信息图](./img/gem_visualization.png)

**细胞分割方法：cellpose V3 模型**

| 质控前 | | 质控后 | |
|--------|--------|--------|--------|
| 细胞数 | {{pre_cell_count}} | 细胞数 | {{post_cell_count}} |
| 中位基因数 | {{pre_median_genes}} | 中位基因数 | {{post_median_genes}} |
| UMI均值 | {{pre_umi}} | UMI均值 | {{post_umi}} |

## 结果展示（开发中）

{{results_content}}

## 交互式结果展示（开发中）


- 细胞类型可视化
- 指定基因表达可视化

*注：完整的交互式结果需要使用专门的可视化工具查看。* 
EOL
        # 创建报告生成脚本
        cat > 04_generate_report.py << 'EOL'
import scanpy as sc
import json
import os
import datetime
import markdown
import sys

def generate_report(chip_id, h5ad_matrix):
    # 读取分析结果
    adata = sc.read_h5ad(h5ad_matrix)
    
    # 准备实验信息数据
    exp_data = {
        "preProcess": {
            "cellCount": str(adata.n_obs),
            "medianGenes": str(int(adata.obs.n_genes.median())),
            "umi": str(int(adata.obs.n_counts.mean()))
        },
        "postProcess": {
            "cellCount": str(len(adata.obs[adata.obs.n_genes > 200])),
            "medianGenes": str(int(adata.obs[adata.obs.n_genes > 200].n_genes.median())),
            "umi": str(int(adata.obs[adata.obs.n_genes > 200].n_counts.mean()))
        }
    }

    # 读取Markdown模板
    with open("report_template.md", "r") as f:
        md_template = f.read()
    
    # 替换模板中的变量
    md_template = md_template.replace("{{patient_name}}", "示例患者")
    md_template = md_template.replace("{{patient_id}}", "示例病理号")
    md_template = md_template.replace("{{sample_date}}", datetime.datetime.now().strftime("%Y-%m-%d"))
    md_template = md_template.replace("{{report_date}}", datetime.datetime.now().strftime("%Y-%m-%d"))
    md_template = md_template.replace("{{clinical_diagnosis}}", "待诊断")
    
    # 实验数据替换
    md_template = md_template.replace("{{pre_cell_count}}", exp_data["preProcess"]["cellCount"])
    md_template = md_template.replace("{{post_cell_count}}", exp_data["postProcess"]["cellCount"])
    md_template = md_template.replace("{{pre_median_genes}}", exp_data["preProcess"]["medianGenes"])
    md_template = md_template.replace("{{post_median_genes}}", exp_data["postProcess"]["medianGenes"])
    md_template = md_template.replace("{{pre_umi}}", exp_data["preProcess"]["umi"])
    md_template = md_template.replace("{{post_umi}}", exp_data["postProcess"]["umi"])
    
    # 添加结果展示占位符处理
    md_template = md_template.replace("{{results_content}}", "结果内容将在后续版本中添加")
    
    # 将Markdown内容转换为HTML
    html_content = markdown.markdown(md_template, extensions=['tables', 'fenced_code'])
    
    # 添加基本的HTML头和CSS样式
    html_doc = f"""<!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>病理诊断模型报告 - {chip_id}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
            h1, h2 {{ color: #333; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            img {{ max-width: 100%; height: auto; }}
        </style>
    </head>
    <body>
    {html_content}
    </body>
    </html>
    """
    
    # 保存最终的HTML文件
    with open("report.html", "w") as f:
        f.write(html_doc)
    
    print(f"报告已生成: report.html")


if __name__ == "__main__":
    # 参数解析
    chip_id = sys.argv[1]
    h5ad_matrix = sys.argv[2]
    
    generate_report(chip_id, h5ad_matrix)
EOL

        # 执行脚本
        python 04_generate_report.py ~{chip_id} ~{h5ad_matrix}
    >>>

    output {
        File report_html = "report.html"
    }

    runtime {
        docker: docker_image
        memory: "16 GB"
        cpu: 4
    }
}