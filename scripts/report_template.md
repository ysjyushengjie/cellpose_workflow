<div style="position: absolute; top: 10px; right: 10px;">
  <img src="./img/BGI_logo.jpg" alt="BGI Logo" style="width: 150px;">
</div>

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