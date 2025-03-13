"""
======================================================
@File   : 04_generate_report.py
@Author : yushengjie
@Date   : 2025.03.13
@Email  : yushengjie@genomics.cn
@Desc   : 从HTML模板生成报告
@Version: v2.0
@Update : 2025.03.14
======================================================
"""

import scanpy as sc
import shutil
import os
import datetime
import sys
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_file_exists(file_path, file_desc):
    """检查文件是否存在"""
    if not os.path.exists(file_path):
        logging.error(f"{file_desc} 文件不存在: {file_path}")
        return False
    return True

def generate_report(chip_id, h5ad_matrix, gem_tiff, html_template):
    """生成HTML报告"""
    # 验证文件存在
    required_files = [
        (h5ad_matrix, "H5AD矩阵文件"),
        (gem_tiff, "GEM TIFF文件"),
        (html_template, "HTML模板文件")
    ]
    
    for file_path, desc in required_files:
        if not check_file_exists(file_path, desc):
            return False
    
    # 创建目录结构
    os.makedirs("img", exist_ok=True)
    os.makedirs("data", exist_ok=True)
    
    # 复制图片文件
    shutil.copy2(gem_tiff, "img/gem_visualization.png")
    
    # 处理BGI logo
    bgi_logo = os.path.join(os.path.dirname(html_template), "BGI_logo.jpg")
    if os.path.exists(bgi_logo):
        shutil.copy2(bgi_logo, "img/BGI_logo.jpg")
    else:
        logging.warning(f"BGI logo文件不存在: {bgi_logo}")
        # 创建一个简单的文本文件作为占位符
        with open("img/BGI_logo.jpg", "w") as f:
            f.write("BGI Logo Placeholder")
    
    # 读取h5ad数据
    try:
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
        
        # 保存实验数据为JSON文件
        with open("data/experimental_data.json", "w", encoding="utf-8") as f:
            json.dump(exp_data, f, ensure_ascii=False, indent=2)
            
    except Exception as e:
        logging.error(f"读取或处理H5AD文件时出错: {str(e)}")
        # 创建默认的实验数据
        exp_data = {
            "preProcess": {
                "cellCount": "N/A",
                "medianGenes": "N/A",
                "umi": "N/A"
            },
            "postProcess": {
                "cellCount": "N/A",
                "medianGenes": "N/A",
                "umi": "N/A"
            }
        }
        with open("data/experimental_data.json", "w", encoding="utf-8") as f:
            json.dump(exp_data, f, ensure_ascii=False, indent=2)
    
    # 生成病人信息数据
    patient_info = {
        "patientName": "示例患者",
        "patientId": "示例病理号",
        "sampleDate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "reportDate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "clinicalDiagnosis": "待诊断"
    }
    
    # 保存病人信息为JSON文件
    with open("data/patient_info.json", "w", encoding="utf-8") as f:
        json.dump(patient_info, f, ensure_ascii=False, indent=2)
    
    # 复制HTML模板到输出文件
    try:
        shutil.copy2(html_template, "report.html")
        logging.info("报告已生成: report.html")
        return True
    except Exception as e:
        logging.error(f"复制HTML模板时出错: {str(e)}")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 5:
        logging.error("参数不足。使用方式: python 04_generate_report.py <chip_id> <h5ad_matrix> <gem_tiff> <html_template>")
        sys.exit(1)
        
    chip_id = sys.argv[1]
    h5ad_matrix = sys.argv[2]
    gem_tiff = sys.argv[3]
    html_template = sys.argv[4]
    
    success = generate_report(chip_id, h5ad_matrix, gem_tiff, html_template)
    sys.exit(0 if success else 1) 