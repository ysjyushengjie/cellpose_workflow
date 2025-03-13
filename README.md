# Cellpose Workflow

基于 Cellpose 的细胞分割和单细胞矩阵生成工作流

## 项目结构

```
.
├── README.md                                    # 项目说明文档
├── cellpose_workflow.wdl                        # WDL工作流主文件
├── cellpose_inputs.json                         # 输入参数配置文件
└── scripts/
    ├── 01_gem2tiff.py                          # GEM转TIFF脚本
    ├── 02_based_cellpose_segment_cell.py       # 细胞分割脚本
    ├── 03_cell_correct_and_generate_matrix.py  # 矩阵生成脚本
    ├── 04_generate_report.py                   # 报告生成脚本
    └── report_template.md                      # 报告Markdown模板
```

## 功能简介

该工作流程实现了以下功能:
1. GEM 文件转换为 TIFF 图像
2. 使用 Cellpose 模型进行细胞分割
3. 生成单细胞矩阵
4. 生成网页可视化报告（基于Markdown模板自动生成HTML）

## 输入参数说明

```
input_tif: 输入的 TIF 图像文件
input_gem: 输入的 GEM 文件
chip_id: 芯片 ID
model_path: Cellpose 模型路径
docker_image: Docker 镜像
scripts_dir: 脚本目录路径
html_template:  Markdown报告模板路径
```

## 输出文件

工作流将生成以下文件:
- `{chip_id}_bin1_beforePS.tiff`: GEM 可视化 TIFF 图像
- `{chip_id}_bin1_beforePS.png`: GEM 可视化 PNG 图像
- `{chip_id}_cellpose_cp_masks.png`: 细胞分割掩码
- `{chip_id}_cellpose.h5ad`: 单细胞矩阵（H5AD格式）
- `{chip_id}.adjusted.cellbin.gef`: 细胞 bin GEF 文件
- `report.html`: 实验报告网页（从Markdown模板生成）

## 报告生成说明

该工作流使用Markdown模板生成HTML报告，而不是直接使用HTML模板。这种方法的优点是:
1. Markdown格式更易于编辑和维护
2. 使用Python标准库可以轻松将Markdown转换为HTML
3. 不需要复杂的JavaScript交互

## 运行方法

1. 修改输入参数:
   编辑 `cellpose_inputs.json` 文件，设置所需参数

2. 运行工作流:
```bash
# 使用 Cromwell 运行
cromwell run cellpose_workflow.wdl -i cellpose_inputs.json

# 或使用 miniwdl
miniwdl run cellpose_workflow.wdl -i cellpose_inputs.json
```

## 系统要求

- Docker
- 16GB 内存
- GPU 支持
- 4 CPU 核心

## Docker 环境

工作流使用预配置的 Docker 镜像，包含所有必要的依赖项:
- Cellpose
- Stereo-seq 相关工具
- Scanpy
- Spateo
- Python Markdown库（用于MD到HTML转换）
