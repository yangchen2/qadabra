import logging
import pandas as pd
from bokeh.plotting import output_file, save
from bokeh.models import DataTable, TableColumn, ColumnDataSource, Div
from bokeh.layouts import column
from bokeh.io import show

# Set up logging
logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

# Set up Bokeh output file
output_file(filename=snakemake.output[0], title="Table")
pval_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

# Process DataFrame
pval_df = pval_df.round(2)
pval_df = pval_df.reset_index()

# Create ColumnDataSource
source = ColumnDataSource(pval_df)

# Define column widths (you can adjust these values as needed)
column_widths = {
    column: 150 for column in pval_df.columns
}

# Create TableColumns with specific widths
columns = [
    TableColumn(field=x, title=x, width=column_widths.get(x, 150))
    for x in pval_df.columns
]

# Create DataTable with sizing_mode
data_table = DataTable(
    source=source,
    columns=columns,
    frozen_columns=1,
    sizing_mode="stretch_width",
    width=1000,
    height=600
)

# Create a title
title = Div(text="<h1>P-value Results Table</h1>", width=1000, height=50)

# Layout the title and the table
layout = column(title, data_table)

# Save the layout to file
save(layout)
logger.info(f"Saved table to {snakemake.output[0]}")