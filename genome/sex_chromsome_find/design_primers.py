#!/01_soft/mambaforge/bin/python
import os
import sys
import getopt
sys.path.append('/01_soft/mambaforge/lib/python3.10/site-packages/')

import pandas as pd
import pysam
import primer3
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="从筛选列表设计 Indel 簇引物")
    parser.add_argument("-i", "--input", required=True, help="输入的 Indel 表格文件 (txt/tsv)")
    parser.add_argument("-g", "--genome", required=True, help="参考基因组 FASTA 文件")
    parser.add_argument("-o", "--out", default="table_designed_primers.tsv", help="输出引物表格")
    
    # 聚类参数
    parser.add_argument("--max_gap", type=int, default=300, help="聚类距离 (bp)，默认300")
    parser.add_argument("--min_net_diff", type=int, default=10, help="最小净长度差异 (bp)")
    return parser.parse_args()

def calculate_delta(ref, alt_str):
    """
    计算 indel 长度差异。
    处理多个 ALT 的情况 (e.g., 'TAC,T')，取差异最大的那个。
    """
    if not isinstance(alt_str, str): return 0
    
    alts = alt_str.split(',')
    best_diff = 0
    
    for alt in alts:
        # 跳过 VCF 中的缺失标记 '*'
        if alt == '*': continue
        
        # 差异 = Alt长度 - Ref长度 (插入为正，缺失为负)
        diff = len(alt) - len(ref)
        
        # 我们取绝对值最大的变异，因为那样跑胶最好看
        if abs(diff) > abs(best_diff):
            best_diff = diff
            
    return best_diff

def design_primer(chrom, start, end, net_diff, genome_fa):
    """调用 Primer3 设计引物"""
    flank = 300 # 往两边各扩 300bp 来找引物
    
    # 提取序列
    fetch_start = max(0, start - flank)
    fetch_end = end + flank
    try:
        seq_template = genome_fa.fetch(chrom, fetch_start, fetch_end)
    except:
        return None

    # 目标区域 (Cluster) 在序列中的位置
    target_start = start - fetch_start
    target_len = end - start
    
    # Primer3 配置
    seq_args = {
        'SEQUENCE_ID': f'{chrom}_{start}',
        'SEQUENCE_TEMPLATE': seq_template,
        'SEQUENCE_TARGET': [target_start, target_len]
    }
    
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[target_len + 50, 800]], # 产物长度
        'PRIMER_NUM_RETURN': 1
    }
    
    try:
        res = primer3.designPrimers(seq_args, global_args)
    except:
        return None
        
    if res['PRIMER_PAIR_NUM_RETURNED'] == 0:
        return None

    ref_size = res['PRIMER_PAIR_0_PRODUCT_SIZE']
    alt_size = ref_size + net_diff
    
    return {
        'Left_Primer': res['PRIMER_LEFT_0_SEQUENCE'],
        'Right_Primer': res['PRIMER_RIGHT_0_SEQUENCE'],
        'Ref_Band': ref_size,
        'Alt_Band': alt_size,
        'Net_Diff': abs(net_diff)
    }

def main():
    args = parse_args()
    
    # 1. 读取表格 (自动处理表头)
    # 假设是用 Tab 分隔的
    try:
        df = pd.read_csv(args.input, sep=None, engine='python')
    except Exception as e:
        print(f"读取表格失败: {e}")
        sys.exit(1)
        
    # 检查列名
    required_cols = ['Chromosome', 'Pos', 'Ref', 'Alt']
    for c in required_cols:
        if c not in df.columns:
            print(f"错误: 表格缺少列 '{c}'。请检查表头。")
            sys.exit(1)
            
    # 2. 准备基因组
    genome = pysam.FastaFile(args.genome)
    
    # 3. 预处理：计算每个位点的 Delta
    df['Delta'] = df.apply(lambda row: calculate_delta(row['Ref'], row['Alt']), axis=1)
    
    # 排序 (聚类前必须排序)
    df = df.sort_values(by=['Chromosome', 'Pos'])
    
    print("正在进行聚类和设计...")
    
    clusters = []
    current_cluster = []
    
    # 4. 聚类逻辑
    for _, row in df.iterrows():
        # 如果 Delta 为 0 (比如 Alt 是 *)，跳过
        if row['Delta'] == 0: continue

        item = {
            'chrom': row['Chromosome'],
            'pos': row['Pos'],
            'ref_len': len(row['Ref']),
            'delta': row['Delta']
        }
        
        if not current_cluster:
            current_cluster.append(item)
            continue
            
        last = current_cluster[-1]
        
        # 判断合并: 同染色体 且 距离 < max_gap
        dist = row['Pos'] - last['pos']
        if row['Chromosome'] == last['chrom'] and dist < args.max_gap:
            current_cluster.append(item)
        else:
            clusters.append(current_cluster)
            current_cluster = [item]
            
    if current_cluster: clusters.append(current_cluster)
    
    print(f"共识别出 {len(clusters)} 个 Indel 候选簇。")
    
    # 5. 设计引物
    results = []
    for cl in clusters:
        # 计算该簇的属性
        chrom = cl[0]['chrom']
        start = cl[0]['pos']
        # End = 最后一个 Indel 的起点 + 它的 Ref 长度
        end = cl[-1]['pos'] + cl[-1]['ref_len']
        
        # 净差异
        net_diff = sum([x['delta'] for x in cl])
        
        # 筛选条件
        if abs(net_diff) < args.min_net_diff:
            continue
            
        res = design_primer(chrom, start, end, net_diff, genome)
        
        if res:
            res['Chromosome'] = chrom
            res['Start'] = start
            res['End'] = end
            res['Indel_Count'] = len(cl)
            results.append(res)
            print(f"成功: {chrom}:{start} (包含 {len(cl)} 个 Indel, 差异 {net_diff}bp)")
            
    # 6. 保存
    if results:
        out_df = pd.DataFrame(results)
        # 整理列顺序
        cols = ['Chromosome', 'Start', 'End', 'Indel_Count', 'Ref_Band', 'Alt_Band', 'Net_Diff', 'Left_Primer', 'Right_Primer']
        out_df = out_df[cols]
        out_df.to_csv(args.out, sep='\t', index=False)
        print(f"\n全部完成！结果已保存至: {args.out}")
        print(out_df.head().to_string(index=False))
    else:
        print("未找到符合条件的引物 (可能是净差异太小)。")

if __name__ == "__main__":
    main()
