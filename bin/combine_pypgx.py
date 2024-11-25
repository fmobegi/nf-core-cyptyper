#!/usr/bin/env python3
import pandas as pd
import glob
import os
import argparse
from openpyxl import load_workbook
from openpyxl.styles import Border, Side, Alignment

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Combine PyPGx results from multiple TSV files."
)
parser.add_argument(
    "-i", "--indir", required=True, help="Input directory containing TSV files"
)
args = parser.parse_args()

# Get all TSV files from the input directory
all_files = glob.glob(os.path.join(args.indir, "*.tsv"))

# Read and combine all TSV files
df_list = []
for filename in all_files:
    df = pd.read_csv(filename, sep="\t", index_col=0)

    # Determine the pipeline based on filename
    if "NGS" in filename:
        pipeline = "NGS"
    elif "TGS" in filename:
        pipeline = "TGS"
    else:
        pipeline = "Unknown"

    # Add the Pipeline column
    df["Pipeline"] = pipeline

    # Reset index to turn the index into a column
    df = df.reset_index()

    # Rename the index column to 'SampleID' if it doesn't have a name
    if df.columns[0] == "index":
        df = df.rename(columns={"index": "SampleID"})

    df_list.append(df)

# Concatenate all dataframes
combined_df = pd.concat(df_list, axis=0)

# Ensure 'SampleID' is the first column
cols = combined_df.columns.tolist()
cols.insert(0, cols.pop(cols.index("SampleID")))
combined_df = combined_df[cols]

# Reorder columns to put Pipeline after Genotype
genotype_index = cols.index("Genotype")
cols.insert(genotype_index + 1, cols.pop(cols.index("Pipeline")))
combined_df = combined_df[cols]

# Write the combined dataframe to a new TSV file
combined_df.to_csv("combined_pypgx_results.tsv", sep="\t", index=False)

# Write the combined dataframe to a new XLSX file
excel_file = "combined_pypgx_results.xlsx"
combined_df.to_excel(excel_file, index=False, engine="openpyxl")

# Load the workbook and worksheet to apply formatting
wb = load_workbook(excel_file)
ws = wb.active

# Define border style
thin_border = Border(
    left=Side(style="thin"),
    right=Side(style="thin"),
    top=Side(style="thin"),
    bottom=Side(style="thin"),
)

# Apply borders to all cells
for row in ws.iter_rows(
    min_row=1, max_row=ws.max_row, min_col=1, max_col=ws.max_column
):
    for cell in row:
        cell.border = thin_border

# Adjust column widths
for column in ws.columns:
    max_length = 0
    column_letter = column[0].column_letter  # Get the column name (e.g., 'A')
    for cell in column:
        try:
            if len(str(cell.value)) > max_length:
                max_length = len(cell.value)
        except:
            pass
    adjusted_width = max_length + 2
    ws.column_dimensions[column_letter].width = adjusted_width

# Save the workbook with formatting applied
wb.save(excel_file)

print(f"Combined results saved to combined_pypgx_results.tsv and {excel_file}")
