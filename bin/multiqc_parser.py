#!/usr/bin/env python3

# Import libraries
import sys
import argparse
import pandas as pd


def parse_multiqc_txt(args):
    
    # Read tab-separated table
    df = pd.read_csv(args.input, sep = "\t")

    # Split column into two columns
    df[["id", "filter_type"]] = df["Sample"].str.rsplit("_", n=1, expand=True)
    # Remove original Sample column
    df.drop("Sample", axis=1, inplace=True)

    # Reshape: long to wide
    df_reads = df.pivot(index="id", columns="filter_type", values="FastQC_mqc_generalstats_fastqc_total_sequences")
    df_length = df.pivot(index="id", columns="filter_type", values="FastQC_mqc_generalstats_fastqc_avg_sequence_length")
    
    # Calculate difference between R1 and R2 raw reads (for paired-end only)
    if "raw.R2" in df_reads.columns:
         df_reads["raw_reads_diff_mln"] = df_reads["raw.R1"] -  df_reads["raw.R2"] 
         df_length["raw_avg_length_diff_bp"] = df_length["raw.R1"] -  df_length["raw.R2"] 
    
    # Round to one decimal point or to the whole number, abs()
    if "raw.R2" in df_reads.columns:
         df_reads["raw_reads_diff_mln"] = abs(df_reads["raw_reads_diff_mln"].astype(int))
    df_reads["raw.R1"] = df_reads["raw.R1"].round(1)
    df_reads["kraken.R1"] = df_reads["kraken.R1"].round(1)
    
    if "raw.R2" in df_reads.columns:
         df_length["raw_avg_length_diff_bp"] = abs(df_length["raw_avg_length_diff_bp"].astype(int))
    df_length["raw.R1"] = df_length["raw.R1"].astype(int)
    df_length["kraken.R1"] = df_length["kraken.R1"].astype(int)
    
    # Reset index
    df_reads = df_reads.reset_index()
    df_length = df_length.reset_index()

    # Remove index name
    df_reads = df_reads.rename_axis(None, axis=1)
    df_length = df_length.rename_axis(None, axis=1)
    
    # Select needed columns
    if "raw.R2" in df_reads.columns:
         df_reads = df_reads.loc[:,["id", "raw_reads_diff_mln", "raw.R1", "kraken.R1"]]
         df_length = df_length.loc[:,["id", "raw_avg_length_diff_bp", "raw.R1", "kraken.R1"]]
    else:
         df_reads = df_reads.loc[:,["id", "raw.R1", "kraken.R1"]]
         df_length = df_length.loc[:,["id", "raw.R1", "kraken.R1"]]
         
    # Rename column names
    df_reads = df_reads.rename(columns = {"raw.R1" : "reads_raw_mln", "kraken.R1" : "reads_final_mln"})
    df_length = df_length.rename(columns = {"raw.R1" : "length_raw", "kraken.R1" : "length_final"})
    
    # Merge df_reads and df_length by a common id column
    df_merged = pd.merge(df_reads, df_length, on="id") 
   
    # Create a clean csv file  
    df_merged.to_csv(args.output, index=False, encoding="utf-8", header=True, sep=",")
   

  

def main():
    parser=argparse.ArgumentParser(description = "Rashape multiqc report table to get a desired format.")
    parser.add_argument("-i", "--input", help="Multiqc report in txt format", dest = "input", type = str, required=True)
    parser.add_argument("-o", '--output', required=True, help="File name with .csv at the end for saving a report.", type = str)
    parser.add_argument("-v", "--verbose", help="verbose", dest="verbose", action='store_true')
    parser.set_defaults(func = parse_multiqc_txt)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()
     
