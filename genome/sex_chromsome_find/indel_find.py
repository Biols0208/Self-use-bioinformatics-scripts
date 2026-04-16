#!/01_soft/mambaforge/bin/python
import os
import sys
import getopt
sys.path.append('/01_soft/mambaforge/lib/python3.10/site-packages/')
import pysam
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="筛选性别特异性 Indel (Presence/Absence Variation)"
    )
    parser.add_argument("-v", "--vcf", required=True, help="输入的 VCF 文件 (建议先只提取 Indel)")
    parser.add_argument("-m", "--males", required=True, help="雄性样本列表")
    parser.add_argument("-f", "--females", required=True, help="雌性样本列表")
    parser.add_argument("-o", "--out", default="sex_specific_indels.txt", help="输出结果文件")
    parser.add_argument("-t", "--tolerance", type=float, default=0.1, 
                        help="容错率 (0.1 表示允许 10%% 的样本不符合规律，默认 0.1)")
    return parser.parse_args()

def get_samples(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def get_gt_status(record, sample_name):
    """
    判断基因型状态:
    0: Reference (0/0)
    1: Variant (0/1, 1/1, 1/2 etc.)
    -1: Missing (./.)
    """
    try:
        gt = record.samples[sample_name]['GT']
        if gt is None or None in gt:
            return -1 # Missing
        
        if gt == (0, 0):
            return 0 # Ref
        else:
            return 1 # Variant (Has Indel)
    except:
        return -1

def main():
    args = parse_args()

    # 1. 读取样本
    raw_males = get_samples(args.males)
    raw_females = get_samples(args.females)
    
    # 2. 打开 VCF
    try:
        vcf_in = pysam.VariantFile(args.vcf)
    except Exception as e:
        print(f"无法打开 VCF 文件: {e}")
        sys.exit(1)

    # 3. 过滤有效样本
    vcf_samples = list(vcf_in.header.samples)
    males = [m for m in raw_males if m in vcf_samples]
    females = [f for f in raw_females if f in vcf_samples]
    
    print(f"分析样本数: 雄性 {len(males)}, 雌性 {len(females)}")
    print(f"容错率: {args.tolerance * 100}%")

    # 4. 打开输出文件
    with open(args.out, 'w') as out_f:
        # 写入表头
        header = ["Chromosome", "Pos", "Ref", "Alt", "Type", "Male_Stats(Var/Tot)", "Female_Stats(Var/Tot)", "Pattern_Score"]
        out_f.write("\t".join(header) + "\n")
        
        count = 0
        passed = 0
        
        # 5. 遍历 VCF
        for record in vcf_in:
            count += 1
            if count % 10000 == 0:
                print(f"已扫描 {count} 个位点...", end="\r")

            # 统计雄性
            m_var = 0
            m_valid = 0
            for m in males:
                st = get_gt_status(record, m)
                if st != -1:
                    m_valid += 1
                    if st == 1: m_var += 1
            
            # 统计雌性
            f_var = 0
            f_valid = 0
            for f in females:
                st = get_gt_status(record, f)
                if st != -1:
                    f_valid += 1
                    if st == 1: f_var += 1

            # 跳过有效样本太少的位点
            if m_valid == 0 or f_valid == 0:
                continue

            # 计算频率
            m_freq = m_var / m_valid
            f_freq = f_var / f_valid
            
            # === 核心筛选逻辑 ===
            match_type = None
            score = 0
            
            # 逻辑 A: 雄性特有 (Male Specific) -> Male High, Female Low
            # 雄性变异率 > (1 - 容错), 雌性变异率 < 容错
            if m_freq >= (1 - args.tolerance) and f_freq <= args.tolerance:
                match_type = "Male_Specific"
                score = m_freq - f_freq # 分数越高越好
            
            # 逻辑 B: 雌性特有 (Female Specific) -> Female High, Male Low
            elif f_freq >= (1 - args.tolerance) and m_freq <= args.tolerance:
                match_type = "Female_Specific"
                score = f_freq - m_freq

            # 如果符合条件，写入文件
            if match_type:
                passed += 1
                
                # 格式化 ALT 列 (可能由多个)
                alt_str = ",".join([str(a) for a in record.alts]) if record.alts else "."
                
                line = [
                    record.chrom,
                    str(record.pos),
                    record.ref,
                    alt_str,
                    match_type,
                    f"{m_var}/{m_valid}({m_freq:.2f})",
                    f"{f_var}/{f_valid}({f_freq:.2f})",
                    f"{score:.2f}"
                ]
                out_f.write("\t".join(line) + "\n")

    print(f"\n扫描完成！")
    print(f"共扫描: {count} 个位点")
    print(f"筛选出: {passed} 个性别特异性 Indel")
    print(f"结果已保存至: {args.out}")

if __name__ == "__main__":
    main()
