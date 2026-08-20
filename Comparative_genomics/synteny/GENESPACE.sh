#!/bin/bash
# =========================================================================
# GENESPACE 完整流程：OrthoFinder + GENESPACE + PAV 矩阵生成
# =========================================================================

# ================= 配置 =================
IN_DIR="/genespace"
SIF_PATH="/01_soft/singularity/genespace.sif"

# 是否跳过已完成的 OrthoFinder（设为 true 则使用已有结果）
SKIP_ORTHOFINDER=true

# ================= 预检 =================
echo "========================================="
echo "GENESPACE Full Pipeline for Sparidae"
echo "========================================="
echo "Start time: $(date)"
echo "Input directory: ${IN_DIR}"
echo "Singularity image: ${SIF_PATH}"
echo ""

# 检查必需目录
if [ ! -d "${IN_DIR}/bed" ] || [ ! -d "${IN_DIR}/peptide" ]; then
    echo "ERROR: 缺少必需目录！"
    echo "Expected:"
    echo "  - ${IN_DIR}/bed"
    echo "  - ${IN_DIR}/peptide"
    exit 1
fi

if [ ! -f "${SIF_PATH}" ]; then
    echo "ERROR:  Singularity 镜像不存在:  ${SIF_PATH}"
    exit 1
fi

# 统计输入文件
BED_COUNT=$(ls ${IN_DIR}/bed/*.bed 2>/dev/null | wc -l)
PEPTIDE_COUNT=$(ls ${IN_DIR}/peptide/*.fa* 2>/dev/null | wc -l)

echo "Input files:"
echo "  BED files: ${BED_COUNT}"
echo "  Peptide files:  ${PEPTIDE_COUNT}"

if [ ${BED_COUNT} -eq 0 ] || [ ${PEPTIDE_COUNT} -eq 0 ]; then
    echo "ERROR: 缺少输入文件！"
    exit 1
fi

echo ""
echo "BED files:"
ls -lh ${IN_DIR}/bed/*.bed
echo ""
echo "Peptide files:"
ls -lh ${IN_DIR}/peptide/*.fa*

# =========================================================================
# 第一步：运行 OrthoFinder
# =========================================================================

if [ -d "${IN_DIR}/orthofinder" ] && [ "${SKIP_ORTHOFINDER}" = "true" ]; then
    echo ""
    echo "========================================="
    echo "跳过 OrthoFinder（使用已有结果）"
    echo "========================================="
    echo "OrthoFinder 目录: ${IN_DIR}/orthofinder"
    
else
    echo ""
    echo "========================================="
    echo "STEP 1: Running OrthoFinder"
    echo "========================================="
    echo ""
    
    # 如果目录存在，先删除
    if [ -d "${IN_DIR}/orthofinder" ]; then
        echo "警告:  orthofinder 目录已存在，将被删除"
        rm -rf ${IN_DIR}/orthofinder
    fi
    
    if [ -d "${IN_DIR}/orthofinder_temp" ]; then
        rm -rf ${IN_DIR}/orthofinder_temp
    fi
    
    echo "运行 OrthoFinder"
    echo "命令:  orthofinder -f peptide -X"
    echo ""
    
    # 运行 OrthoFinder
    /usr/bin/singularity exec --cleanenv \
        -B ${IN_DIR}:/input \
        ${SIF_PATH} \
        orthofinder \
            -f /input/peptide \
            -t 48 \
            -a 48 \
            -X \
            -o /input/orthofinder_temp
    
    OF_EXIT=$?
    
    echo ""
    echo "OrthoFinder finished with exit code: ${OF_EXIT}"
    
    if [ ${OF_EXIT} -ne 0 ]; then
        echo "ERROR: OrthoFinder 失败！"
        exit ${OF_EXIT}
    fi
    
    # 移动结果到标准位置
    RESULT_DIR=$(find ${IN_DIR}/orthofinder_temp -maxdepth 1 -type d -name "Results_*" | head -1)
    
    if [ -z "$RESULT_DIR" ]; then
        echo "ERROR: 找不到 OrthoFinder 结果目录！"
        exit 1
    fi
    
    mv ${IN_DIR}/orthofinder_temp ${IN_DIR}/orthofinder
    echo ""
    echo "✓ OrthoFinder 完成！"
fi

# =========================================================================
# 第二步：生成 R 脚本
# =========================================================================

echo ""
echo "========================================="
echo "STEP 2: Preparing GENESPACE"
echo "========================================="

R_SCRIPT="${IN_DIR}/run_genespace_full.R"
echo "生成 R 脚本: ${R_SCRIPT}"

cat > ${R_SCRIPT} <<'RSCRIPT'
library(GENESPACE)

# ========== 日志设置 ==========
cat("=========================================\n")
cat("GENESPACE Multi-Species Analysis\n")
cat("=========================================\n")
cat(paste("Start time:", Sys.time(), "\n"))
cat(paste("Working directory:", getwd(), "\n"))
cat(paste("R version:", R.version.string, "\n"))
cat(paste("GENESPACE version:", packageVersion("GENESPACE"), "\n\n"))

# ========== 初始化 GENESPACE ==========
cat("========== Initializing GENESPACE ==========\n")
gpar <- init_genespace(
    wd = "/input",
    path2mcscanx = "/opt/MCScanX",
    rawOrthofinderDir = "/input/orthofinder"
)

cat("\nDetected genomes:\n")
print(gpar$genomeIDs)

# ========== 运行 GENESPACE ==========
cat("\n========== Running GENESPACE Pipeline ==========\n")
cat("注意:  将使用已有的 OrthoFinder 结果\n\n")

out <- run_genespace(gpar, overwrite = FALSE)

cat("\n========== Pipeline Complete ==========\n")

# ========== 生成 PAV 矩阵 ==========
cat("\n========== Generating PAV Matrices ==========\n")

all_genomes <- gpar$genomeIDs
cat(paste("Total genomes:", length(all_genomes), "\n"))

out_dir <- "/input/PAV_Matrices"
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

# 压平函数
flatten_df <- function(df) {
    df[] <- lapply(df, function(col) {
        if (is.list(col)) {
            sapply(col, function(x) {
                if (is.null(x) || length(x) == 0 || all(is.na(x))) {
                    return(NA_character_)
                }
                paste(unique(x), collapse = ";")
            })
        } else {
            return(col)
        }
    })
    return(df)
}

success_count <- 0
fail_count <- 0

for (i in seq_along(all_genomes)) {
    ref_id <- all_genomes[i]
    cat(paste0("\n[", i, "/", length(all_genomes), "] Processing:  ", ref_id, "\n"))
    
    tryCatch({
        pg_matrix <- query_pangenes(gpar, refGenome = ref_id, transform = TRUE)
        
        if (is.null(pg_matrix) || nrow(pg_matrix) == 0) {
            warning(paste("No data for", ref_id))
            fail_count <- fail_count + 1
            next
        }
        
        df_out <- flatten_df(as.data.frame(pg_matrix))
        out_file <- file.path(out_dir, paste0("PAV_matrix_Ref_", ref_id, ".csv"))
        write.csv(df_out, out_file, row.names = FALSE)
        
        cat(paste("  ✓ Saved:", basename(out_file), "\n"))
        cat(paste("    Size:", nrow(df_out), "rows x", ncol(df_out), "cols\n"))
        
        success_count <- success_count + 1
        
    }, error = function(e) {
        cat(paste("  ✗ ERROR:", e$message, "\n"))
        fail_count <- fail_count + 1
    })
}

# ========== 汇总统计 ==========
cat("\n=========================================\n")
cat("Final Summary\n")
cat("=========================================\n")
cat(paste("PAV matrices generated:", success_count, "/", length(all_genomes), "\n"))
cat(paste("Failed:", fail_count, "\n"))

if (success_count > 0) {
    summary_file <- file.path(out_dir, "PAV_summary.txt")
    sink(summary_file)
    cat("GENESPACE Analysis Summary\n")
    cat("==========================\n\n")
    cat(paste("Date:", Sys.time(), "\n"))
    cat(paste("Genomes:", length(all_genomes), "\n"))
    cat(paste("Genome IDs:", paste(all_genomes, collapse = ", "), "\n\n"))
    
    if (!is.null(out$pangenes)) {
        cat(paste("Total orthogroups:", length(unique(out$pangenes$og)), "\n\n"))
    }
    
    cat("Gene counts:\n")
    for (genome in all_genomes) {
        count <- sum(gpar$synteny$bed$genome == genome)
        cat(paste("  ", genome, ":", count, "\n"))
    }
    sink()
    
    cat(paste("\nSummary saved to:", summary_file, "\n"))
}

cat("\n=========================================\n")
cat(paste("End time:", Sys.time(), "\n"))
cat("All done!\n")
cat("=========================================\n")
RSCRIPT
# =========================================================================
# 第三步：运行 GENESPACE
# =========================================================================

echo ""
echo "========================================="
echo "STEP 3: Running GENESPACE"
echo "========================================="
echo ""
echo "开始运行 GENESPACE..."
echo ""

# 运行 GENESPACE
/usr/bin/singularity exec --cleanenv \
    -B ${IN_DIR}:/input \
    ${SIF_PATH} \
    Rscript /input/run_genespace_full.R

GS_EXIT=$?

# =========================================================================
# 最终报告
# =========================================================================

echo ""
echo "========================================="
echo "Pipeline Finished"
echo "========================================="
echo "End time: $(date)"
echo "Exit code: ${GS_EXIT}"
echo ""

if [ ${GS_EXIT} -eq 0 ]; then
    echo "✓ 流程成功完成！"
    echo ""
    
    # 显示结果摘要
    echo "========================================="
    echo "Results Summary"
    echo "========================================="
    
    if [ -d "${IN_DIR}/orthofinder" ]; then
        echo ""
        echo "OrthoFinder 结果:"
        if [ -f "${IN_DIR}/orthofinder/Orthogroups/Orthogroups.tsv" ]; then
            OG_COUNT=$(tail -n +2 ${IN_DIR}/orthofinder/Orthogroups/Orthogroups.tsv | wc -l)
            echo "  Orthogroups:  ${OG_COUNT}"
        fi
    fi
    
    if [ -d "${IN_DIR}/PAV_Matrices" ]; then
        echo ""
        echo "PAV 矩阵文件:"
        ls -lh ${IN_DIR}/PAV_Matrices/*.csv 2>/dev/null
    fi
    
    if [ -d "${IN_DIR}/results" ]; then
        echo ""
        echo "主要结果 (${IN_DIR}/results/):"
        ls -lh ${IN_DIR}/results/ 2>/dev/null | head -10
    fi
    
    if [ -d "${IN_DIR}/riparian" ]; then
        echo ""
        echo "可视化图表 (${IN_DIR}/riparian/):"
        ls -lh ${IN_DIR}/riparian/*.pdf 2>/dev/null | head -5
    fi
    
    echo ""
    echo "========================================="
    echo "完成！所有结果保存在:"
    echo "${IN_DIR}"
    echo "========================================="
    
else
    echo "✗ 流程失败！"
    echo ""
    echo "请检查错误日志:"
    echo "  genespace_full_${SLURM_JOB_ID}.err"
    echo ""
    echo "可能的问题:"
    echo "  1. 检查 GENESPACE 日志输出"
    echo "  2. 确认 MCScanX 路径正确"
    echo "  3. 检查内存是否足够"
fi

exit ${GS_EXIT}
