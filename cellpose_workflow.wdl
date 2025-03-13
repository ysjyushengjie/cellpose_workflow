version 1.0

# Cellpose 工作流：用于细胞分割和单细胞矩阵生成
# Cellpose Workflow: For cell segmentation and single-cell matrix generation
workflow cellpose_workflow {
    input {
        # 输入的 TIF 图像文件 | Input TIF image file
        File input_tif
        # 输入的 GEM 文件 | Input GEM file
        File input_gem
        # 芯片 ID | Chip ID
        String chip_id
        # Cellpose 模型路径 | Cellpose model path
        String model_path
        # Docker 镜像地址 | Docker image address
        String docker_image
        # 脚本目录路径 | Scripts directory path
        String scripts_dir
        # HTML模板路径 | HTML template path
        String html_template
        
        # 运行环境配置参数 | Runtime configuration parameters
        Int memory_gb = 16
        Int cpu_count = 4
        Int gpu_count = 1
        String gpu_type = "L4"
    }

    # GEM转换为TIFF | Convert GEM to TIFF
    call gem2tiff {
        input:
            input_gem = input_gem,
            chip_id = chip_id,
            docker_image = docker_image,
            scripts_dir = scripts_dir,
            memory_gb = memory_gb,
            cpu_count = cpu_count,
            gpu_count = gpu_count,
            gpu_type = gpu_type
    }

    # 调用细胞分割任务 | Call cell segmentation task
    call segment_cell {
        input:
            input_tif = input_tif,
            model_path = model_path,
            chip_id = chip_id,
            docker_image = docker_image,
            scripts_dir = scripts_dir,
            memory_gb = memory_gb,
            cpu_count = cpu_count,
            gpu_count = gpu_count,
            gpu_type = gpu_type
    }

    # 调用矩阵生成任务 | Call matrix generation task
    call generate_matrix {
        input:
            input_gem = input_gem,
            chip_id = chip_id,
            mask_path = segment_cell.mask_output,
            docker_image = docker_image,
            scripts_dir = scripts_dir,
            memory_gb = memory_gb,
            cpu_count = cpu_count,
            gpu_count = gpu_count,
            gpu_type = gpu_type
    }

    # 生成实验报告 | Generate report
    call generate_report {
        input:
            chip_id = chip_id,
            h5ad_matrix = generate_matrix.h5ad_output,
            gem_tiff = gem2tiff.tiff_output,
            docker_image = docker_image,
            html_template = html_template,
            scripts_dir = scripts_dir,
            memory_gb = memory_gb,
            cpu_count = cpu_count
    }

    output {
        File gem_tiff = gem2tiff.tiff_output
        File gem_png = gem2tiff.png_output
        # 细胞掩码文件 | Cell mask file
        File mask_file = segment_cell.mask_output
        # H5AD 格式的单细胞矩阵 | Single-cell matrix in H5AD format
        File h5ad_matrix = generate_matrix.h5ad_output
        # 细胞 bin GEF 文件 | Cell bin GEF file
        File cell_bin_gef = generate_matrix.cell_bin_gef
        # 网页报告 | HTML report
        File html_report = generate_report.report_html
    }
}

# 新增 GEM 转 TIFF 任务
task gem2tiff {
    input {
        File input_gem
        String chip_id
        String docker_image
        String scripts_dir
        Int memory_gb
        Int cpu_count
        Int gpu_count
        String gpu_type
    }

    command <<<
        python ~{scripts_dir}/01_gem2tiff.py ~{input_gem} ~{chip_id}
    >>>

    output {
        File tiff_output = "${chip_id}_bin1_beforePS.tiff"
        File png_output = "${chip_id}_bin1_beforePS.png"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu_count
        gpu: gpu_count
        gpu_type: gpu_type
    }
}

# 细胞分割任务：使用 Cellpose 模型进行细胞分割
# Cell Segmentation Task: Using Cellpose model for cell segmentation
task segment_cell {
    input {
        # 输入的 TIF 图像文件 | Input TIF image file
        File input_tif
        # Cellpose 模型路径 | Cellpose model path
        String model_path
        # 芯片 ID | Chip ID
        String chip_id
        # Docker 镜像地址 | Docker image address
        String docker_image
        # 脚本目录路径 | Scripts directory path
        String scripts_dir
        Int memory_gb
        Int cpu_count
        Int gpu_count
        String gpu_type
    }

    # 执行命令：调用 Python 脚本进行细胞分割
    # Execute command: Call Python script for cell segmentation
    command <<<
        python ~{scripts_dir}/02_based_cellpose_segment_cell.py ~{input_tif} ~{model_path} ~{chip_id}
    >>>

    # 任务输出：细胞掩码文件 | Task output: Cell mask file
    output {
        File mask_output = "${chip_id}_cellpose_cp_masks.png"
    }

    # 运行时配置：使用 GPU 加速 | Runtime configuration: Using GPU acceleration
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu_count
        gpu: gpu_count
        gpu_type: gpu_type
    }
}

# 细胞矩阵生成任务：基于分割结果生成单细胞矩阵
# Cell Matrix Generation Task: Generate single-cell matrix based on segmentation results
task generate_matrix {
    input {
        # 输入的 GEM 文件 | Input GEM file
        File input_gem
        # 芯片 ID | Chip ID
        String chip_id
        # 细胞掩码文件路径 | Cell mask file path
        File mask_path
        # Docker 镜像地址 | Docker image address
        String docker_image
        # 脚本目录路径 | Scripts directory path
        String scripts_dir
        Int memory_gb
        Int cpu_count
        Int gpu_count
        String gpu_type
    }

    # 执行命令：调用 Python 脚本生成单细胞矩阵
    # Execute command: Call Python script to generate single-cell matrix
    command <<<
        python ~{scripts_dir}/03_cell_correct_and_generate_matrix.py ~{input_gem} ~{chip_id} ~{mask_path}
    >>>

    # 任务输出：H5AD 矩阵和 GEF 文件 | Task output: H5AD matrix and GEF file
    output {
        File h5ad_output = "${chip_id}_cellpose.h5ad"
        File cell_bin_gef = "${chip_id}.adjusted.cellbin.gef"
    }

    # 运行时配置 | Runtime configuration
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu_count
        gpu: gpu_count
        gpu_type: gpu_type
    }
}

task generate_report {
    input {
        String chip_id
        File h5ad_matrix
        File gem_tiff
        String docker_image
        String html_template
        String scripts_dir
        Int memory_gb
        Int cpu_count
    }

    command <<<
        python ~{scripts_dir}/04_generate_report.py \
            ~{chip_id} \
            ~{h5ad_matrix} \
            ~{gem_tiff} \
            ~{scripts_dir}/index.html
    >>>

    output {
        File report_html = "report.html"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu_count
    }
}
