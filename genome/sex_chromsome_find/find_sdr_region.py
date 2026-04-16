#!/01_soft/mambaforge/bin/python
import os
import sys
import getopt
sys.path.append('/01_soft/mambaforge/lib/python3.10/site-packages/')

import pandas as pd
import numpy as np
import sys

# ================= 配置区域 =================
INPUT_CSV = "sex_heterozygosity_diff.csv"
OUTPUT_CANDIDATES = "candidate_SDR_regions.csv"
OUTPUT_BED = "candidate_SDR.bed"  # 方便导入 IGV 查看

# 阈值设置
# 1. 差异阈值：Het_Diff 绝对值大于多少才算有效信号？(建议 0.6 - 0.8)
DIFF_THRESHOLD = 0.8

# 2. 合并距离 (bp)：两个显著点相距多远以内算作同一个区域？
# 鱼类 SDR 可能很小，也可能很大.
MAX_GAP = 10000 

# 3. 最小 SNP 数：一个区域至少要有几个支持点才算数？(防止单个测序错误)
MIN_SNPS = 3
# ===========================================

def find_sdr_candidates(file_path):
    print(f"正在读取数据: {file_path} ...")
    try:
        df = pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        sys.exit(1)

    # 1. 筛选显著位点
    # 关心的是 |Diff| > 阈值 的点
    sig_df = df[df['Het_Diff'].abs() >= DIFF_THRESHOLD].copy()
    
    if sig_df.empty:
        print("未找到任何显著差异位点，请尝试降低 DIFF_THRESHOLD。")
        return

    # 确保按染色体和位置排序
    sig_df = sig_df.sort_values(by=['CHR', 'POS'])
    
    print(f"筛选出 {len(sig_df)} 个显著差异 SNP (阈值 > {DIFF_THRESHOLD})")
    print("正在进行区域合并 (Clustering)...")

    regions = []
    
    # 初始化第一个簇
    if len(sig_df) > 0:
        # 获取第一行数据
        first_row = sig_df.iloc[0]
        current_region = {
            'CHR': first_row['CHR'],
            'Start': first_row['POS'],
            'End': first_row['POS'],
            'SNPs': [first_row['Het_Diff']], # 记录差异值用于算平均
            'Count': 1
        }

    # 遍历剩余的行
    for i in range(1, len(sig_df)):
        row = sig_df.iloc[i]
        
        # 判断是否属于同一个区域：
        # 1. 染色体相同
        # 2. 位置距离 < MAX_GAP
        if (row['CHR'] == current_region['CHR']) and \
           (row['POS'] - current_region['End'] <= MAX_GAP):
            
            # 更新当前区域
            current_region['End'] = row['POS']
            current_region['SNPs'].append(row['Het_Diff'])
            current_region['Count'] += 1
            
        else:
            # 结束上一个区域，保存
            regions.append(current_region)
            
            # 开启新区域
            current_region = {
                'CHR': row['CHR'],
                'Start': row['POS'],
                'End': row['POS'],
                'SNPs': [row['Het_Diff']],
                'Count': 1
            }
    
    # 不要忘了保存最后一个区域
    if len(sig_df) > 0:
        regions.append(current_region)

    # 2. 整理结果并推断性别系统
    results = []
    for r in regions:
        # 过滤掉 SNP 数量太少的区域
        if r['Count'] < MIN_SNPS:
            continue
            
        avg_diff = np.mean(r['SNPs'])
        length = r['End'] - r['Start'] + 1
        
        # 推断类型
        # Diff < 0 (负值): 雄性杂合度高 -> 同型XY系统的SDR，或ZW系统的Z
        # Diff > 0 (正值): 雌性杂合度高 -> 异型XY系统的X，或同型ZW系统的SDR
        inference = "Unknown"
        if avg_diff < -0.5:
            inference = "Likely XY SDR (Homomorphic) / ZW-Z"
        elif avg_diff > 0.5:
            inference = "Likely XY-X (Heteromorphic) / ZW SDR"

        results.append({
            'Chromosome': r['CHR'],
            'Start': r['Start'],
            'End': r['End'],
            'Length_bp': length,
            'SNP_Count': r['Count'],
            'Mean_Het_Diff': round(avg_diff, 4),
            'Inference': inference
        })

    # 转为 DataFrame 输出
    res_df = pd.DataFrame(results)
    
    if res_df.empty:
        print("没有区域满足最小 SNP 数量要求。")
        return

    # 按 SNP 数量和平均差异强度排序
    res_df['Abs_Mean_Diff'] = res_df['Mean_Het_Diff'].abs()
    res_df = res_df.sort_values(by=['SNP_Count', 'Abs_Mean_Diff'], ascending=False)
    res_df = res_df.drop(columns=['Abs_Mean_Diff']) # 移除辅助列

    # 保存 CSV
    res_df.to_csv(OUTPUT_CANDIDATES, index=False, sep="\t")
    print(f"\n成功识别 {len(res_df)} 个候选区域！已保存至: {OUTPUT_CANDIDATES}")
    
    # 保存 BED 文件 (用于 IGV)
    # BED 格式: chrom start end name score
    with open(OUTPUT_BED, 'w') as f:
        for _, row in res_df.iterrows():
            # BED start is 0-based, output is 1-based, strictly speaking minus 1, but for approximation it's fine.
            # Name format: Diff_Count
            name = f"Diff:{row['Mean_Het_Diff']}|cnt:{row['SNP_Count']}"
            
            # 设置颜色: 红色为正差异(XY-X)，蓝色为负差异(XY-SDR)
            color = "255,0,0" if row['Mean_Het_Diff'] > 0 else "0,0,255"
            
            # 标准 BED9 格式便于 IGV 展示颜色
            # chrom start end name score strand thickStart thickEnd itemRgb
            line = f"{row['Chromosome']}\t{row['Start']-1}\t{row['End']}\t{name}\t{row['SNP_Count']}\t.\t{row['Start']-1}\t{row['End']}\t{color}\n"
            f.write(line)
            
    print(f"BED 文件已保存至: {OUTPUT_BED} (可拖入 IGV 查看)")
    print("-" * 50)
    print(res_df.head(10).to_string(index=False))

if __name__ == "__main__":
    find_sdr_candidates(INPUT_CSV)
