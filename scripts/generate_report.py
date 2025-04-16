"""
======================================================
@File   : generate_report.py
@Author : yushengjie
@Date   : 2025.03.13
@Email  : yushengjie@genomics.cn
@Desc   : Generate HTML report from h5ad matrix and image files
@mirror : spateo_scanpy
@Version: v4.0
@Update : 2025.04.15
======================================================

Update record:
    1. Added cell ratio image to HTML report
    2. Added distribution image to HTML report
    3. Added cnv result image to HTML report
    4. Added log file to record the process
    5. Added error handling
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
    try:
        if not os.path.exists(file_path):
            logging.error(f"error: {file_desc} does not exist at {file_path}")
            return False
        logging.info(f"{file_desc} found at {file_path}")
        return True
    except Exception as e:
        logging.error(f"error checking file {file_path}: {str(e)}")
        return False

def generate_experimental_data(h5ad_matrix):
    """
    generate experimental data from h5ad matrix

    Parameters:
        h5ad_matrix (str): H5AD matrix file path
    Returns:
        dict: experimental data
    """
    try:
        logging.info(f"Loading h5ad matrix from {h5ad_matrix}")
        adata = sc.read_h5ad(h5ad_matrix)
        logging.info(f"Successfully loaded h5ad matrix with {adata.n_obs} cells and {adata.n_vars} genes")
        
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

def generate_report(html_template, logo_path, chip_id, h5ad_matrix, gem_png, cell_ratio_image, distribution_image, cnv_result_image):
    """
    generate html report from h5ad matrix and image files

    Parameters:
        html_template (str): HTML template file path
        logo_path (str): Logo file path

        chip_id (str): chip id
        h5ad_matrix (str): H5AD matrix file path
        gem_png (str): GEM PNG file path

        cell_ratio_image (str): cell ratio image path
        distribution_image (str): cell distribution image path
        cnv_result_image (str): CNV result image path

    Return:
        bool: True if report generation is successful, False otherwise
    """
    logging.info(f"Starting report generation for chip {chip_id}")
    logging.info(f"Input files: h5ad={h5ad_matrix}, gem={gem_png}, template={html_template}")

    # check if the required files exist
    required_files = [
        (html_template, "HTML template file path"),
        (logo_path, "BGI logo file path"),

        (h5ad_matrix, "H5AD matrix file path"),
        (gem_png, "GEM PNG file path"),

        (cell_ratio_image, "cell ratio image path"),
        (distribution_image, "cell distribution image path"),
        (cnv_result_image, "CNV result image path"),
    ]
    logging.info("<==========================loading required files========================>")
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
        logging.info("<======================generating report=======================>")
        with open(html_template, "r", encoding="utf-8") as template_file:
            if not template_file:
                logging.error(f"error: HTML template file {html_template} is empty")
                return False
            html_content = template_file.read()
            logging.info(f"HTML template loaded successfully from {html_template}")

        # generate embedded data
        logging.info("generating embedded data")
        embedded_patient_info = f"<script id='embedded-patient-info'>const patientInfo = {json.dumps(patient_info, ensure_ascii=False)};</script>"
        embedded_exp_data = f"<script id='embedded-exp-data'>const experimentalData = {json.dumps(exp_data, ensure_ascii=False)};</script>"

        # replace placeholders or insert into HTML header
        html_content = html_content.replace("{{LOGO_PATH}}", logo_path)
        html_content = html_content.replace("{{GEM_IMAGE_PATH}}", gem_png)
        html_content = html_content.replace("{{CELL_RATIO_IMAGE}}", cell_ratio_image)
        html_content = html_content.replace("{{DISTRIBUTION_IMAGE}}", distribution_image)
        html_content = html_content.replace("{{CNV_RESULT_IMAGE}}", cnv_result_image)
        html_content = html_content.replace("<!-- EMBEDDED_DATA_PLACEHOLDER -->", 
                                            f"{embedded_patient_info}\n{embedded_exp_data}")

        with open("report.html", "w", encoding="utf-8") as output_file:
            if not output_file:
                logging.error(f"error: output file {output_file} is empty")
                return False
            output_file.write(html_content)

        logging.info("report generation successful")
        return True
    except Exception as e:
        logging.error(f"error generating report: {str(e)}")
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--html", required=True, help="HTML")
    parser.add_argument("--logo", default="img/BGI_logo.jpg", help="LOGO")

    parser.add_argument("--chip", required=True, help="chipId")
    parser.add_argument("--h5ad", required=True, help="H5AD")
    parser.add_argument("--gem", required=True, help="PNG")

    parser.add_argument("--cell_ratio", default="img/cell_ratio.png", help="cell ratio image path")
    parser.add_argument("--distribution", default="img/distribution.png", help="cell distribution image path")
    parser.add_argument("--cnv_result", default="img/cnv_result.png", help="CNV result image path")

    args = parser.parse_args()

    success = generate_report(
        html_template = args.html,
        logo_path=args.logo,

        chip_id = args.chip,
        h5ad_matrix = args.h5ad,

        gem_png = args.gem,
        cell_ratio_image=args.cell_ratio,
        distribution_image=args.distribution,
        cnv_result_image=args.cnv_result
    )
