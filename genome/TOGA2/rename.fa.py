#!/ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/miniconda3/bin/python

"""
特点:
1. 基于BioPython的高效FASTA解析
2. 如未找到映射则保持原样不变
3. 匹配成功时保留原ID: >new old other_information
4. 支持序列中的特殊字符如'-'
5. 支持制表符或空格分隔的映射文件
6. 使用精确匹配（区分大小写）进行ID对比
"""

import os
import sys
import getopt
sys.path.append('/ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/miniconda3/lib/python3.7/site-packages/')
import argparse
from argparse import RawTextHelpFormatter
import re

# 简单用法信息
SIMPLE_USAGE = """
用法: rename_fasta.py [--ids old_new.txt] input.fasta [--out rename.fasta]

示例:
  rename_fasta.py --ids new_names.txt input.fasta > renamed.fasta
  rename_fasta.py --ids names_map.txt input.fasta --out output.fasta
"""

# 函数定义
def msg(*args, **kwargs):
    """向stderr输出信息"""
    print(*args, file=sys.stderr, **kwargs)

def err(*args, **kwargs):
    """向stderr输出错误信息并退出"""
    msg(*args, **kwargs)
    sys.exit(1)

# 检查文件是否存在
def check_file(f):
    return os.path.isfile(f)

# 检查文件是否为FASTA格式
def quick_check_fasta(f):
    if not os.path.isfile(f) or os.path.getsize(f) < 1:
        return False
    with open(f, 'r') as fasta:
        if fasta.readline()[0] != '>':  # 检查头部是否以">"开始
            return False
    return True

# 加载ID映射表，支持制表符或空格分隔
def load_id_mapping(mapping_file):
    """灵活加载ID映射，支持制表符或空格分隔"""
    id_map = {}
    line_num = 0
    
    try:
        with open(mapping_file, 'r', encoding='utf-8-sig') as f:  # 处理可能的BOM
            for line in f:
                line_num += 1
                line = line.strip()
                
                # 跳过空行和注释行
                if not line or line.startswith('#'):
                    continue
                
                # 尝试多种分隔符 - 首先尝试制表符，然后尝试空格
                if '\t' in line:
                    parts = line.split('\t')
                else:
                    # 使用正则表达式分割一个或多个空格
                    parts = re.split(r'\s+', line, 1)
                
                if len(parts) < 2:
                    msg(f"警告: 第{line_num}行格式不正确 '{line}'，需要至少两列")
                    continue
                
                old_id = parts[0].strip()
                new_id = parts[1].strip()
                
                # 验证ID有效性
                if not old_id:
                    msg(f"警告: 第{line_num}行原始ID为空")
                    continue
                    
                if not new_id:
                    msg(f"警告: 第{line_num}行新ID为空，将使用原始ID")
                    new_id = old_id
                
                # 存储映射关系 - 保持原始大小写，区分大小写
                id_map[old_id] = new_id
    
    except Exception as e:
        err(f"解析ID映射文件失败: {e}")
    
    if not id_map:
        err(f"未从'{mapping_file}'中加载到有效的ID映射")
        
    return id_map

def process_fasta(fasta_file, id_map, output_file=None, verbose=False):
    """处理FASTA文件 - 仅在需要时导入BioPython"""
    # 导入BioPython模块 - 仅在实际处理文件时加载
    try:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from io import StringIO
    except ImportError:
        err("无法导入BioPython。请安装： pip install biopython")
    
    # 处理FASTA文件
#    msg(f"处理FASTA文件...")
    newseqs = []
    stats = {"total": 0, "renamed": 0, "unchanged": 0}
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        stats["total"] += 1
        original_id = record.id
        
        # 原始完整描述（去除ID部分）
        full_desc = record.description[len(original_id):].strip()
        
        # 查找映射
        if original_id in id_map:
            new_id = id_map[original_id]
            
            if verbose:
                msg(f"重命名: '{original_id}' -> '{new_id}'")
            
            # 创建新序列记录
            newseqs.append(SeqRecord(record.seq, id=new_id, description=f"{original_id}{' ' + full_desc if full_desc else ''}"))
            stats["renamed"] += 1
        else:
            # 未找到映射，保持原样
            if verbose:
                msg(f"保持不变: '{original_id}'")
                
            newseqs.append(SeqRecord(record.seq, id=original_id, description=full_desc))
            stats["unchanged"] += 1
    
    # 写入文件或输出到标准输出
    if output_file:
        SeqIO.write(newseqs, output_file, 'fasta')
        msg(f"处理完成: 已保存到 {output_file}")
    else:
        seqFILE = StringIO()
        SeqIO.write(newseqs, seqFILE, 'fasta')
        output = seqFILE.getvalue().rstrip()
        print(output)
    
    # 打印统计信息
#    msg(f"处理完成: 总共 {stats['total']} 条序列, 重命名 {stats['renamed']} 条, 保持不变 {stats['unchanged']} 条")
    
    return stats

def main():
    # 检查是否仅需要显示帮助
    if len(sys.argv) <= 1 or '-h' in sys.argv or '--help' in sys.argv:
        print(SIMPLE_USAGE)
        return 0
    
    # 命令行参数解析
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description='FASTA序列标识符重命名工具 - 性能优化版\n' +
                    '1. 如未找到映射则保持原样不变\n' +
                    '2. 匹配成功时保留原ID: >new old other_information\n' +
                    '3. 支持制表符或空格分隔的映射文件\n' +
                    '4. 使用精确匹配（区分大小写）进行ID对比\n',
        usage='\n  %(prog)s [--ids new_names.txt] FASTA > new.fasta')
    parser.add_argument('fasta', metavar='FASTA', nargs=1, help='原始FASTA文件')
    parser.add_argument('--ids', metavar='FILE', required=True, nargs=1, 
                        help='指定包含[旧名称][制表符或空格][新名称]的文件')
    parser.add_argument('--out', metavar='FILE', nargs=1, help='指定输出文件（默认为标准输出）')
    parser.add_argument('--verbose', action='store_true', help='显示详细处理信息')
    parser.add_argument('--version', action='version', version='%(prog)s v0.3.0')
    args = parser.parse_args()
    
    # 打印运行信息
#    msg(f"-- 运行时间: {CURRENT_TIME} --")
#    msg(f"-- 运行用户: {CURRENT_USER} --")
    
    # 是否显示详细信息
    verbose = args.verbose
    
    # 检查输入/输出文件
    if not check_file(args.fasta[0]):
        err(f'ERROR: 找不到文件"{args.fasta[0]}"，请检查指定目录中是否存在该文件。')
    if not quick_check_fasta(args.fasta[0]):
        err(f'ERROR: 请检查"{args.fasta[0]}"是否为FASTA格式。')
    if not check_file(args.ids[0]):
        err(f'ERROR: 找不到文件"{args.ids[0]}"，请检查指定目录中是否存在该文件。')
    if args.out and check_file(args.out[0]):
        err(f'ERROR: 输出文件"{args.out[0]}"已存在。')
    
    # 加载ID映射
#    msg(f"加载ID映射文件...")
    id_map = load_id_mapping(args.ids[0])
#    msg(f"已加载 {len(id_map)} 个ID映射关系")
    
    # 处理FASTA文件 - 输出文件名提取
    output_file = args.out[0] if args.out else None
    
    # 调用处理函数
    process_fasta(args.fasta[0], id_map, output_file, verbose)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
