import sys
import statistics

# --- Configuration for Fluctuation Detection ---
FLUCTUATION_STD_DEV_THRESHOLD = 1.5
# --- End Configuration ---

def calculate_stats(numbers_str_list):
    """
    Calculates statistics for a list of number strings.
    Returns: A list of strings for all calculated statistics.
    """
    # Initialize all metrics to default "N/A" or "0" for error or empty cases
    # Order: Zeros, Sum, Mean, Median, Variance, StdDev, Fluctuations,
    #        Min, Max, CV(%), Present_Count, Present_%, Single_Copy_Count, Single_Copy_%
    default_return_values = ["0", "0", "N/A", "N/A", "N/A", "N/A", "0", # Base + Fluctuations
                             "N/A", "N/A", "N/A", "0", "N/A", "0", "N/A"] # New stats

    if not numbers_str_list: # No numbers after colon
        return default_return_values

    try:
        numbers = [int(n) for n in numbers_str_list]
    except ValueError: # Invalid number format
        error_return = ["INVALID_NUM_FORMAT"] * len(default_return_values)
        return error_return

    if not numbers: # Should be caught by 'if not numbers_str_list', but as a safeguard
        return default_return_values

    # --- Existing Calculations ---
    zeros_count = numbers.count(0)
    total_sum = sum(numbers)
    
    mean_val_float = None
    stdev_val_float = None
    
    mean_val_str = "N/A"
    median_val_str = "N/A"
    variance_val_str = "N/A" 
    stdev_val_str = "N/A"
    num_sharp_fluctuations_str = "N/A"

    try:
        mean_val_float = statistics.mean(numbers)
        mean_val_str = f"{mean_val_float:.2f}"
    except statistics.StatisticsError: pass 

    try:
        median_val = statistics.median(numbers)
        median_val_str = f"{median_val:.2f}" if isinstance(median_val, float) else str(median_val) 
    except statistics.StatisticsError: pass

    if len(numbers) >= 2:
        try:
            variance_val_str = f"{statistics.variance(numbers):.2f}"
            stdev_val_float = statistics.stdev(numbers)
            stdev_val_str = f"{stdev_val_float:.2f}"
        except statistics.StatisticsError:
            variance_val_str = "CALC_ERROR"; stdev_val_str = "CALC_ERROR"; stdev_val_float = None 
    elif len(numbers) == 1:
        variance_val_str = "0.00"; stdev_val_str = "0.00"; stdev_val_float = 0.0

    if mean_val_float is not None and stdev_val_float is not None:
        if stdev_val_float == 0:
            num_sharp_fluctuations_str = "0"
        else:
            fluctuations_count = sum(1 for num in numbers if abs(num - mean_val_float) > FLUCTUATION_STD_DEV_THRESHOLD * stdev_val_float)
            num_sharp_fluctuations_str = str(fluctuations_count)
    elif len(numbers) == 1 and mean_val_float is not None:
        num_sharp_fluctuations_str = "0"
    else:
        num_sharp_fluctuations_str = "N/A"
        if variance_val_str == "CALC_ERROR": num_sharp_fluctuations_str = "N/A"

    # --- New Calculations ---
    # 1. Range
    min_count_str = str(min(numbers)) if numbers else "N/A"
    max_count_str = str(max(numbers)) if numbers else "N/A"

    # 2. Coefficient of Variation (CV)
    cv_str = "N/A"
    if mean_val_float is not None and stdev_val_float is not None:
        if mean_val_float == 0:
            cv_str = "N/A" # Or "0.00%" if stdev is also 0, but N/A for undefined case
        else:
            cv_str = f"{(stdev_val_float / mean_val_float) * 100:.2f}"
    
    # 3. Presence Breadth
    total_species = len(numbers)
    species_present_count_val = sum(1 for x in numbers if x > 0)
    species_present_count_str = str(species_present_count_val)
    species_present_percent_str = f"{(species_present_count_val / total_species) * 100:.2f}" if total_species > 0 else "N/A"

    # 4. Single-Copy Statistics
    single_copy_species_count_val = sum(1 for x in numbers if x == 1)
    single_copy_species_count_str = str(single_copy_species_count_val)
    single_copy_species_percent_str = f"{(single_copy_species_count_val / total_species) * 100:.2f}" if total_species > 0 else "N/A"

    return [
        str(zeros_count), str(total_sum), mean_val_str, median_val_str, variance_val_str, stdev_val_str, num_sharp_fluctuations_str,
        min_count_str, max_count_str, cv_str, 
        species_present_count_str, species_present_percent_str,
        single_copy_species_count_str, single_copy_species_percent_str
    ]

def process_file(input_filepath, output_filepath):
    processed_lines_output = []
    first_line_processed = False
    stats_headers = [
        "Zeros", "Sum", "Mean", "Median", "Variance", "StdDev", f"Fluctuations(>{FLUCTUATION_STD_DEV_THRESHOLD}std)",
        "Min_Count", "Max_Count", "CV(%)", 
        "Species_Present_Count", "Species_Present_Percent(%)",
        "Single_Copy_Species_Count", "Single_Copy_Species_Percent(%)"
    ]

    try:
        with open(input_filepath, 'r', encoding='utf-8') as infile:
            for line_num, raw_line in enumerate(infile, 1):
                line_content = raw_line.strip()
                output_line = line_content 

                if not line_content: 
                    processed_lines_output.append("")
                    if not first_line_processed: first_line_processed = True
                    continue

                if not first_line_processed:
                    if not line_content.startswith("OG"): 
                        output_line = f"{line_content} {' '.join(stats_headers)}"
                    first_line_processed = True

                if line_content.startswith("OG") and ':' in line_content:
                    parts = line_content.split(':', 1)
                    numbers_data_str = parts[1]
                    numbers_str_list = numbers_data_str.strip().split() 
                                                                  
                    stats_results = calculate_stats(numbers_str_list)

                    if line_content.startswith("OG"): 
                        if "INVALID_NUM_FORMAT" in stats_results:
                            output_line = f"{line_content} # ERROR: Invalid number format in data at line {line_num}"
                        elif "CALC_ERROR" in stats_results : # CALC_ERROR might be in variance/stdev
                            # Check if any result is CALC_ERROR, signifying an issue
                            is_calc_error = any(res == "CALC_ERROR" for res in stats_results)
                            if is_calc_error:
                                output_line = f"{line_content} # ERROR: Calculation error for some stats at line {line_num}"
                            else: # If no specific CALC_ERROR, join normally
                                output_line = f"{line_content} {' '.join(stats_results)}"
                        else:
                            if len(stats_results) == len(stats_headers):
                                output_line = f"{line_content} {' '.join(stats_results)}"
                            else: 
                                output_line = f"{line_content} # ERROR: Mismatch in stats results length ({len(stats_results)} vs {len(stats_headers)})"
                
                elif line_content.startswith("OG"): 
                    if line_num > 1 or (line_num == 1 and line_content.startswith("OG")):
                         output_line = f"{line_content} # ERROR: Malformed OG line at line {line_num}"
                
                processed_lines_output.append(output_line)
        
        with open(output_filepath, 'w', encoding='utf-8') as outfile:
            for out_line in processed_lines_output:
                outfile.write(out_line + "\n")
        
        print(f"Processing complete. Results written to {output_filepath}")
        print(f"Fluctuation detection used a threshold of {FLUCTUATION_STD_DEV_THRESHOLD} standard deviations.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filepath}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <input_file> <output_file>")
        print("Example: python your_script_name.py data.txt results.txt")
        sys.exit(1)

    # It's good practice to replace 'your_script_name.py' in usage with actual script name
    # For example, by using sys.argv[0] or a fixed name if you rename the file.
    # For now, the placeholder is fine.

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    process_file(input_filename, output_filename)
