"""
======================================================
@File   : 04_generate_report.py
@Author : yushengjie
@Date   : 2025.03.13
@Email  : yushengjie@genomics.cn
@Desc   : generate html report from h5ad matrix and gem png
@mirror : spateo_scanpy
@Version: v3.0
@Update : 2025.03.21
======================================================

Update record:
    1. add embedded data to html header
    2. add logo to html header
"""

import scanpy as sc
import numpy as np
import os
import datetime
import argparse
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_file_exists(file_path, file_desc):
    """
    check if the file exists
    
    Parameters: 
        file_path (str): file path
        file_desc (str): file description
    Returns:
        bool: True if the file exists, False otherwise
    """
    if not os.path.exists(file_path):
        logging.error(f"error: {file_desc} does not exist")
        return False
    return True

def generate_experimental_data(h5ad_matrix):
    """
    generate experimental data from h5ad matrix

    Parameters:
        h5ad_matrix (str): H5AD matrix file path
    Returns:
        dict: experimental data
    """
    try:
        adata = sc.read_h5ad(h5ad_matrix)
        
        # check if the raw data exists
        if adata.raw is not None:
            n_genes_field = "n_genes" if "n_genes" in adata.raw.obs.columns else None
            n_counts_field = "n_counts" if "n_counts" in adata.raw.obs.columns else None
        else:
            # calculate gene and UMI
            n_genes_field = "n_genes"
            n_counts_field = "n_counts"
            adata.obs[n_genes_field] = np.sum(adata.X > 0, axis=1).A1
            adata.obs[n_counts_field] = np.sum(adata.X, axis=1).A1

        return {
            "preProcess": {
                "cellCount": str(adata.n_obs),
                "medianGenes": str(int(adata.obs[n_genes_field].median())) if n_genes_field else "N/A",
                "umi": str(int(adata.obs[n_counts_field].mean())) if n_counts_field else "N/A"
            },
            "postProcess": {
                "cellCount": str(len(adata.obs[adata.obs[n_genes_field] > 200])) if n_genes_field else "N/A",
                "medianGenes": str(int(adata.obs[adata.obs[n_genes_field] > 200][n_genes_field].median())) if n_genes_field else "N/A",
                "umi": str(int(adata.obs[adata.obs[n_genes_field] > 200][n_counts_field].mean())) if n_counts_field else "N/A"
            }
        }
        
    except Exception as e:
        logging.error(f"error generating experimental data: {str(e)}")
        return {
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

def generate_report(chip_id, h5ad_matrix, gem_png, html_template, logo_path):
    """
    generate html report from h5ad matrix and gem png

    Parameters:
        chip_id (str): chip id
        h5ad_matrix (str): H5AD matrix file path
        gem_png (str): GEM PNG file path
        html_template (str): HTML template file path
        logo_path (str): Logo file path

    Return:
        bool: True if report generation is successful, False otherwise
    """

    # check if the required files exist
    required_files = [
        (h5ad_matrix, "H5AD matrix file path"),
        (gem_png, "GEM PNG file path"),
        (html_template, "HTML template file path"),
        (logo_path, "BGI logo file path")
    ]
    
    for file_path, desc in required_files:
        if not check_file_exists(file_path, desc):
            return False

    # generate experimental data
    exp_data = generate_experimental_data(h5ad_matrix)
    
    # generate patient info
    patient_info = {
        "patientName": "示例患者",
        "patientId": "示例病理号",
        "sampleDate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "reportDate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "clinicalDiagnosis": "待诊断"
    }
    
    # generate html report
    try:
        with open(html_template, "r", encoding="utf-8") as template_file:
            html_content = template_file.read()

        # generate embedded data
        embedded_patient_info = f"<script id='embedded-patient-info'>const patientInfo = {json.dumps(patient_info, ensure_ascii=False)};</script>"
        embedded_exp_data = f"<script id='embedded-exp-data'>const experimentalData = {json.dumps(exp_data, ensure_ascii=False)};</script>"

        # replace placeholders or insert into HTML header
        html_content = html_content.replace("{{LOGO_PATH}}", logo_path)
        html_content = html_content.replace("{{GEM_IMAGE_PATH}}", gem_png)
        html_content = html_content.replace("<!-- EMBEDDED_DATA_PLACEHOLDER -->", 
                                            f"{embedded_patient_info}\n{embedded_exp_data}")

        with open("report.html", "w", encoding="utf-8") as output_file:
            output_file.write(html_content)

        logging.info("report generation successful")
        return True
    except Exception as e:
        logging.error(f"error generating report: {str(e)}")
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--chip", required=True, help="chipId")
    parser.add_argument("--h5ad", required=True, help="H5AD")
    parser.add_argument("--gem", required=True, help="PNG")
    parser.add_argument("--html", required=True, help="HTML")
    parser.add_argument("--logo", default="img/BGI_logo.jpg", help="LOGO")

    args = parser.parse_args()

    success = generate_report(
        chip_id = args.chip,
        h5ad_matrix = args.h5ad,
        gem_png = args.gem,
        html_template = args.html,
        logo_path=args.logo
    )
    # args = {
    #     "chip" : "Y00723M3",
    #     "h5ad" : "Y00723M3_cellpose_NC.h5ad",
    #     "gem" : "Y00723M3_bin1_beforePS.png",
    #     "html" : "index.html",
    #     "logo" : "img/BGI_logo.jpg"
    # }
    # success = generate_report(
    #     chip_id = args['chip'],  # 改为方括号访问
    #     h5ad_matrix = args['h5ad'],
    #     gem_png = args['gem'],
    #     html_template = args['html'],
    #     logo_path=args['logo']
    # )
