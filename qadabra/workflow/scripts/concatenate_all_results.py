import logging
import pandas as pd
import io

# Setup logger
logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(f"[%(asctime)s - {snakemake.rule}] :: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

# Function to remove the second row if it contains strings
def remove_string_row(filepath, logger):
    logger.info(f"Checking for string content in the second row of file: {filepath}")
    with open(filepath, "r") as file:
        lines = file.readlines()

    if len(lines) > 1:
        second_line = lines[1].strip().split("\t")
        if any(not cell.replace(".", "", 1).replace("-", "", 1).isdigit() for cell in second_line):
            logger.info("Second row contains strings; removing it.")
            lines.pop(1)  # Remove the second line

    # Save the cleaned content back to a temporary DataFrame
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)

# Load and clean differentials
logger.info("Loading and cleaning differentials...")
concatenated_diff_file = remove_string_row(snakemake.input["concatenated_differentials"], logger)

for col in concatenated_diff_file.columns:
    new_col_name = col + " differentials"
    concatenated_diff_file.rename(columns={col: new_col_name}, inplace=True)

# Load p-values
logger.info("Loading p-values...")
concatenated_pvalue_file = pd.read_csv(snakemake.input["concatenated_pvalues"], sep="\t", index_col=0)

for col in concatenated_pvalue_file.columns:
    new_col_name = col + " significant features"
    concatenated_pvalue_file.rename(columns={col: new_col_name}, inplace=True)

# Concatenate results
logger.info("Concatenating results...")
qadabra_all_result = pd.concat([concatenated_diff_file, concatenated_pvalue_file], axis=1)

# Process significant features
logger.info("Processing significant features...")
for col in qadabra_all_result.columns:
    if col.endswith(" significant features"):
        qadabra_all_result[col] = qadabra_all_result[col].astype(object)  # Ensure compatibility for strings

for index, row in qadabra_all_result.iterrows():
    count = 0
    for col in qadabra_all_result.columns:
        if col.endswith(" significant features"):
            if row[col] < 0.05:
                count += 1
                qadabra_all_result.at[index, col] = "p<0.05"
            else:
                qadabra_all_result.at[index, col] = "ns"
    qadabra_all_result.at[index, "# of tools p > 0.05"] = count

# Save results
logger.info(f"Saving final results to {snakemake.output[0]}...")
qadabra_all_result.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info("Processing complete!")
