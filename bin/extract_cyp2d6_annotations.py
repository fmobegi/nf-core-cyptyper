#!/usr/bin/env python3

import os
import re
import pandas as pd
import argparse
from collections import defaultdict
from openpyxl import load_workbook
from openpyxl.styles import Border, Side, Alignment

# # Function to extract CYP2D6 information and barcode from an HTML file
# def extract_cyp2d6_info(file_path):
#     gene = "CYP2D6"
#     file_name = os.path.basename(file_path)
#     try:
#         # sample_name_match = re.search(r"^(.*?)_barcode(\d+)_", file_name)
#         sample_name_match = re.search(r"^(.*?)_barcode(\d+)[_.]", file_name)
#         if sample_name_match:
#             sample_name = sample_name_match.group(1)
#             barcode = sample_name_match.group(2)
#         else:
#             sample_name = "Unknown"
#             barcode = "Unknown"

#         population_match = re.search(r"_population_(\w+)", file_path)
#         population = population_match.group(1) if population_match else "Unknown"

#         with open(file_path, "r") as file:
#             content = file.read()
#             match = re.search(
#                 rf"<p><font .*?><b>Gene</b>: {gene}.*?<b>Diplotype</b>: (.*?)&nbsp;.*?<b>Phenotype</b>: (.*?)</font>",
#                 content,
#                 re.DOTALL,
#             )
#             if match:
#                 diplotype = match.group(1).strip()
#                 phenotype = match.group(2).strip()
#                 # Remove the word "Metabolizer" from phenotype
#                 phenotype = re.sub(r"\bMetabolizer\b", "", phenotype).strip()
#                 sample_name_with_barcode = f"{sample_name}_barcode{barcode}"
#                 return (
#                     sample_name_with_barcode,
#                     population,
#                     f"{diplotype} | {phenotype}",
#                 )
#     except Exception as e:
#         print(f"Error extracting information from {file_path}: {e}")
#     return f"Unknown_barcode{barcode}", "Unknown", None


# # Parse command-line arguments
# parser = argparse.ArgumentParser(
#     description="Extract CYP2D6 information from HTML files."
# )
# parser.add_argument(
#     "-i", "--input_directory", required=True, help="Directory containing HTML files."
# )
# args = parser.parse_args()

# input_directory = args.input_directory

# try:
#     input_files = [
#         os.path.join(input_directory, f)
#         for f in os.listdir(input_directory)
#         if f.endswith(".html")
#     ]
# except Exception as e:
#     print(f"Error listing files in directory {input_directory}: {e}")
#     input_files = []

# cyp2d6_data = defaultdict(lambda: defaultdict(str))

# for file_path in input_files:
#     try:
#         sample_name, population, diplotype_phenotype = extract_cyp2d6_info(file_path)
#         if diplotype_phenotype:
#             cyp2d6_data[sample_name][population] = diplotype_phenotype
#     except Exception as e:
#         print(f"Error processing file {file_path}: {e}")

# populations = [
#     "Oceania",
#     "AfricanAmerican",
#     "American",
#     "EasternAsia",
#     "European",
#     "Latin",
#     "NearEastern",
#     "SouthASIA",
#     "SubSaharanAfrica",
# ]

# data = []
# for sample_name, population_data in cyp2d6_data.items():
#     row = [sample_name]
#     for pop in populations:
#         row.append(population_data.get(pop, ""))
#     data.append(row)

# columns = ["SampleID"] + populations
# df = pd.DataFrame(data, columns=columns)

# # Drop columns without data
# df = df.dropna(axis=1, how="all")

# # Add citation row
# citation_row = [
#     "Table citation: For each population, results are Diplotype | Phenotype"
# ] + [""] * len(populations)
# df.loc[len(df)] = citation_row

# output_excel_file = "CYP2D6_annotations.xlsx"

# # Write DataFrame to Excel
# try:
#     df.to_excel(output_excel_file, index=False)
# except Exception as e:
#     print(f"Error writing DataFrame to Excel file {output_excel_file}: {e}")

# # Load the workbook and worksheet to apply formatting
# try:
#     wb = load_workbook(output_excel_file)
#     ws = wb.active

#     # Define border style
#     thin_border = Border(
#         left=Side(style="thin"),
#         right=Side(style="thin"),
#         top=Side(style="thin"),
#         bottom=Side(style="thin"),
#     )

#     # Apply borders to all cells
#     for row in ws.iter_rows(
#         min_row=1, max_row=ws.max_row, min_col=1, max_col=ws.max_column
#     ):
#         for cell in row:
#             cell.border = thin_border

#     # Adjust column widths
#     for column in ws.columns:
#         max_length = 0
#         column_letter = column[0].column_letter  # Get the column name (e.g., 'A')
#         for cell in column:
#             try:
#                 if len(str(cell.value)) > max_length:
#                     max_length = len(cell.value)
#             except Exception as e:
#                 print(f"Error calculating cell length: {e}")
#                 pass
#         adjusted_width = max_length + 2
#         ws.column_dimensions[column_letter].width = adjusted_width

#     # Add citation formatting
#     citation_cell = ws.cell(row=ws.max_row, column=1)
#     citation_cell.alignment = Alignment(wrap_text=True, horizontal="left")

#     # Save the workbook with formatting applied
#     wb.save(output_excel_file)
# except Exception as e:
#     print(f"Error formatting or saving Excel file {output_excel_file}: {e}")

# print(f"CYP2D6 population data extracted, formatted, and saved to {output_excel_file}")


## ~~~~~~~~~~~~~~~~~~ ATTEMPT TO ONLY USE POPULATIONS GENOTYPED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Function to extract CYP2D6 information and barcode from an HTML file
def extract_cyp2d6_info(file_path):
    gene = "CYP2D6"
    file_name = os.path.basename(file_path)
    try:
        sample_name_match = re.search(r"^(.*?)_barcode(\d+)[_.]", file_name)
        if sample_name_match:
            sample_name = sample_name_match.group(1)
            barcode = sample_name_match.group(2)
        else:
            sample_name = "Unknown"
            barcode = "Unknown"

        population_match = re.search(r"_population_(\w+)", file_path)
        population = population_match.group(1) if population_match else "Unknown"

        with open(file_path, "r") as file:
            content = file.read()
            match = re.search(
                rf"<p><font .*?><b>Gene</b>: {gene}.*?<b>Diplotype</b>: (.*?)&nbsp;.*?<b>Phenotype</b>: (.*?)</font>",
                content,
                re.DOTALL,
            )
            if match:
                diplotype = match.group(1).strip()
                phenotype = match.group(2).strip()
                phenotype = re.sub(r"\bMetabolizer\b", "", phenotype).strip()
                sample_name_with_barcode = f"{sample_name}_barcode{barcode}"
                return (
                    sample_name_with_barcode,
                    population,
                    f"{diplotype} | {phenotype}",
                )
    except Exception as e:
        print(f"Error extracting information from {file_path}: {e}")
    return f"Unknown_barcode{barcode}", "Unknown", None


# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Extract CYP2D6 information from HTML files."
)
parser.add_argument(
    "-i", "--input_directory", required=True, help="Directory containing HTML files."
)
args = parser.parse_args()

input_directory = args.input_directory

try:
    input_files = [
        os.path.join(input_directory, f)
        for f in os.listdir(input_directory)
        if f.endswith(".html")
    ]
except Exception as e:
    print(f"Error listing files in directory {input_directory}: {e}")
    input_files = []

cyp2d6_data = defaultdict(lambda: defaultdict(str))

for file_path in input_files:
    try:
        sample_name, population, diplotype_phenotype = extract_cyp2d6_info(file_path)
        if diplotype_phenotype:
            cyp2d6_data[sample_name][population] = diplotype_phenotype
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

# Get the list of populations present in the data
populations_present = sorted(
    set(pop for sample_data in cyp2d6_data.values() for pop in sample_data.keys())
)

data = []
for sample_name, population_data in cyp2d6_data.items():
    row = [sample_name]
    for pop in populations_present:
        row.append(population_data.get(pop, ""))
    data.append(row)

columns = ["SampleID"] + populations_present
df = pd.DataFrame(data, columns=columns)

# Add citation row
citation_row = [
    "Table citation: For each population, results are Diplotype | Phenotype"
] + [""] * len(populations_present)
df.loc[len(df)] = citation_row

output_excel_file = "CYP2D6_annotations.xlsx"

# Write DataFrame to Excel
try:
    df.to_excel(output_excel_file, index=False)
except Exception as e:
    print(f"Error writing DataFrame to Excel file {output_excel_file}: {e}")

# Load the workbook and worksheet to apply formatting
try:
    wb = load_workbook(output_excel_file)
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
            except Exception as e:
                print(f"Error calculating cell length: {e}")
                pass
        adjusted_width = max_length + 2
        ws.column_dimensions[column_letter].width = adjusted_width

    # Add citation formatting
    citation_cell = ws.cell(row=ws.max_row, column=1)
    citation_cell.alignment = Alignment(wrap_text=True, horizontal="left")

    # Save the workbook with formatting applied
    wb.save(output_excel_file)
except Exception as e:
    print(f"Error formatting or saving Excel file {output_excel_file}: {e}")

print(f"CYP2D6 population data extracted, formatted, and saved to {output_excel_file}")
