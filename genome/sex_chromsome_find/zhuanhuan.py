#! /01_software/miniconda3/bin/python3
import sys
import csv

def convert_scientific_to_decimal(input_file, output_file, columns):
    with open(input_file, 'r') as f_in:
        reader = csv.reader(f_in, delimiter='\t')
        header = next(reader)
        has_header = csv.Sniffer().has_header('\t'.join(header))
        if has_header:
            next(reader)
        with open(output_file, 'w') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            writer.writerow(header)
            for row in reader:
                for col_idx in columns:
                    try:
                        row[col_idx - 1] = format(float(row[col_idx - 1]), '.10f')
                    except ValueError:
                        pass
                writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py input_file output_file column_index1 column_index2 ...")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    columns = [int(idx) for idx in sys.argv[3:]]
    
    convert_scientific_to_decimal(input_file, output_file, columns)

