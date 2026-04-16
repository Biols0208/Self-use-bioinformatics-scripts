#!/01_soft/mambaforge/bin/python
import os
import sys
import getopt
sys.path.append('/01_soft/mambaforge/lib/python3.10/site-packages/')

import pandas as pd
import matplotlib.pyplot as plt

# ================= 配置区域 =================
MALE_FILE = "males_hardy.hwe"
FEMALE_FILE = "females_hardy.hwe"
OUTPUT_CSV = "sex_heterozygosity_diff.csv"
OUTPUT_IMG = "sex_heterozygosity_plot.png"
# ===========================================

def parse_vcftools_hwe(file_path, label):
    """
    读取 .hwe 文件并计算观测杂合度
    """
    print(f"正在读取文件: {file_path} ...")
    try:
        df = pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        sys.exit(1)

    # 核心函数：解析 "10/5/12" 这种格式
    def calculate_obs_het(obs_str):
        try:
            # obs_str 示例: "10/5/12" (HomRef/Het/HomAlt)
            counts = [int(x) for x in obs_str.split('/')]
            total = sum(counts)
            if total == 0:
                return 0.0
            # 杂合度 = 杂合子数量 / 总样本数
            return counts[1] / total
        except:
            return 0.0

    # 应用函数
    # VCFTools 的列名通常是 'OBS(HOM1/HET/HOM2)'
    col_name = 'OBS(HOM1/HET/HOM2)'
    if col_name not in df.columns:
        print(f"错误: 在 {file_path} 中找不到列 {col_name}")
        # 尝试查找类似的列名（防止版本差异）
        for c in df.columns:
            if 'OBS' in c:
                col_name = c
                break
    
    df[f'Het_{label}'] = df[col_name].apply(calculate_obs_het)
    
    # 只保留需要的列
    return df[['CHR', 'POS', f'Het_{label}']]

# 1. 读取并处理数据
males = parse_vcftools_hwe(MALE_FILE, "Male")
females = parse_vcftools_hwe(FEMALE_FILE, "Female")

# 2. 合并数据 (基于 染色体CHR 和 位置POS)
print("正在合并雌雄数据...")
merged = pd.merge(males, females, on=['CHR', 'POS'], how='inner')

# 3. 计算差异 (Diff = Female - Male)
# 如果是 XY 系统: X染色体上 Female(1.0) - Male(0.0) > 0 (正值)
# 如果是 ZW 系统: Z染色体上 Female(0.0) - Male(1.0) < 0 (负值)
merged['Het_Diff'] = merged['Het_Female'] - merged['Het_Male']

# 4. 排序并保存
merged = merged.sort_values(by=['CHR', 'POS'])
merged.to_csv(OUTPUT_CSV, index=False, sep='\t')
print(f"结果已保存至: {OUTPUT_CSV}")

# 5. 简单统计与筛选
# 筛选差异显著的位点 (比如差异绝对值大于 0.5)
candidates = merged[merged['Het_Diff'].abs() > 0.5]
print("\nTop 候选差异位点预览:")
print(candidates.head())

# ================= 简单的绘图 (Manhattan style) =================
print("\n正在绘图...")

# 为了画图，如果染色体不是数字，需要把它们转为分类索引
merged['CHR_Code'] = pd.factorize(merged['CHR'])[0]

plt.figure(figsize=(15, 6))

# 画背景点 (灰色)
plt.scatter(merged.index, merged['Het_Diff'], c='lightgray', s=2, alpha=0.5, label='Insignificant')

# 画显著点 (红色)
# 阈值设为 0.5 (即两性杂合度差异超过 50%)
significant = merged[merged['Het_Diff'].abs() > 0.5]
plt.scatter(significant.index, significant['Het_Diff'], c='red', s=10, alpha=0.8, label='Significant Diff > 0.5')

plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('Genomic Position (Cumulative Index)')
plt.ylabel('Heterozygosity Difference (Female - Male)')
plt.title('Sex-Specific Heterozygosity Difference')
plt.legend()

# 添加 XY/ZW 提示
plt.text(0, -0.8, "Diff > 0: Likely XY System (XX vs XY)", color='blue', fontsize=10)
plt.text(0, 0.8, "Diff < 0: Likely ZW System (ZW vs ZZ)", color='green', fontsize=10)

plt.tight_layout()
plt.savefig(OUTPUT_IMG, dpi=300)
print(f"图片已保存至: {OUTPUT_IMG}")
