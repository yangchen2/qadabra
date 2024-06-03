import logging
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

plt.style.use(snakemake.config["stylesheet"])

logger.info(f"Reading in concatenated coefficients table")
# Read in coefficients table
concat_coef = pd.read_csv(snakemake.input[0], sep="\t")

logger.info(f"Reading in concatenated pvalues table")
# Read in p-values table
concat_pvalue = pd.read_csv(snakemake.input[1], sep="\t")



logger.info(f"Summarizing results")
# Create column that counts number of tools with significant p-value
concat_pvalue['num_sig'] = concat_pvalue.apply(lambda row: sum(row[1:] < 0.01), axis=1)

# Merging concat_coef with concat_pvalue to add the 'num_sig' column based on the 'feature_id' column
concat_coef = pd.merge(concat_coef, concat_pvalue[['feature_id', 'num_sig']], on='feature_id', how='left')

# Average coefficient across tools and create avg_coef
numeric_columns = concat_coef.columns[1:-1]
concat_coef['avg_coef'] = concat_coef[numeric_columns].mean(axis=1)
df_sorted = concat_coef.sort_values(by='avg_coef', ascending=False)

# Standard deviation and standard error of mean columns
df_sorted['stand_dev'] = df_sorted.iloc[:, 1:-2].std(axis=1)
df_sorted['SEM'] = df_sorted['stand_dev'] / np.sqrt(8)


logger.info(f"Plotting figure")
num_top_features = min(20, len(df_sorted[df_sorted['avg_coef'] > 0]))
num_bottom_features = min(20, len(df_sorted[df_sorted['avg_coef'] < 0]))

top_features = df_sorted.nlargest(num_top_features, 'avg_coef')
bottom_features = df_sorted.nsmallest(num_bottom_features, 'avg_coef')
combined = pd.concat([top_features, bottom_features])
combined_sorted = combined.sort_values(by='avg_coef', ascending=False)  # Sort to have positives on top



logger.info(f"Reading in BIOM table to get abundance info")
# Read in BIOM table
biom_table = load_table(snakemake.input[2])
df = pd.DataFrame(biom_table.to_dataframe())

# Get row sum
df['row_sum'] = df.sum(axis=1)
# Logarithmic normalization
df['Log_Abundance'] = np.log10(df['row_sum'] + 1)  # Adding 1 to avoid log(0)

# Min-Max normalization on the log-transformed data
df['Normalized_Log_Abundance'] = (df['Log_Abundance'] - df['Log_Abundance'].min()) / (df['Log_Abundance'].max() - df['Log_Abundance'].min())

# Convert to dense array before applying the minimum threshold
df['Normalized_Log_Abundance'] = df['Normalized_Log_Abundance'].astype(np.float64)
df['Normalized_Log_Abundance'] = df['Normalized_Log_Abundance'].apply(lambda x: max(x, 0.05))


# Reset index
df = df.reset_index().rename(columns={'index': 'feature_id'})

# Merging both dataframes
combined_sorted = combined_sorted.merge(df[['feature_id', 'Normalized_Log_Abundance']], on='feature_id', how='left')



logger.info(f"Plotting figure")
# def adjust_colors(value, max_positive, max_negative, alpha):
#     if value > 0:
#         return (1, 0, 0, alpha)  # Red with alpha
#     else:
#         return (0, 0, 1, alpha)  # Blue with alpha

def adjust_colors(value, max_positive, max_negative, alpha):
    if value > 0:
        return (1, 0.5, 0.31, alpha)  # Coral with alpha
    else:
        return (0, 0.5, 0.5, alpha)  # Teal with alpha

# Assuming combined_sorted is your DataFrame
max_positive = combined_sorted['avg_coef'][combined_sorted['avg_coef'] > 0].max()
max_negative = -combined_sorted['avg_coef'][combined_sorted['avg_coef'] < 0].min()
alpha_values = combined_sorted['Normalized_Log_Abundance']

# Adjust colors with alpha based on Normalized_Log_Abundance
colors_adjusted = combined_sorted.apply(
    lambda row: adjust_colors(row['avg_coef'], max_positive, max_negative, row['Normalized_Log_Abundance']),
    axis=1
)

# Create a colormap with 8 discrete colors from white to black
colors_8_reversed = np.linspace(1, 0, 8)  # Generate 8 grayscale values from white to black
cmap_8_colors_reversed = ListedColormap([(c, c, c, 1.0) for c in colors_8_reversed], name='custom_grey_8_reversed')

# Create figure and GridSpec layout with adjusted width ratios
fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 0.5])

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Remove top and right spines on the left subplot
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Ensure 'xerr' is a 1D array
xerr = combined_sorted['SEM'].values if combined_sorted['SEM'].shape == combined_sorted['avg_coef'].shape else None

if xerr is None:
    raise ValueError("The shape of 'xerr' does not match the shape of 'avg_coef'. Please check the data.")

# Left Subplot: Horizontal barplot
y_positions = np.arange(len(combined_sorted))
ax1.barh(y_positions, combined_sorted['avg_coef'], color=colors_adjusted.tolist(), xerr=combined_sorted['SEM'].loc[combined_sorted.index])

ax1.set_yticks(y_positions)
ax1.set_yticklabels(combined_sorted['feature_id'], fontsize=8)
ax1.set_xlabel('Average Coefficient')
# ax1.set_title('Summary Plot', fontsize=13)

# ax1.set_title('Summary Plot')

ax1.invert_yaxis()  # Positive values on top

# Right Subplot: Heatmap with 8 reversed discrete colors, making cells square
img = ax2.imshow(combined_sorted['num_sig'].values.reshape(-1, 1), cmap=cmap_8_colors_reversed, aspect='auto', vmin=0, vmax=7)
ax2.set_xticks([])
ax2.set_yticks(y_positions)
ax2.set_yticklabels([])  # Y-tick labels are removed as they're aligned now

# Synchronize y-axis limits across subplots to ensure alignment
ax1.set_ylim(ax2.get_ylim())

# Add colorbar
cbar = fig.colorbar(img, ax=ax2, orientation='vertical', ticks=np.arange(8), label='Number of Tools Adjusted P-value < 0.01')


# Create colormaps from white to blue and white to red
cmap_white_to_blue = LinearSegmentedColormap.from_list('white_to_blue', [(1, 1, 1), (0, 0, 1)], N=100)
cmap_white_to_red = LinearSegmentedColormap.from_list('white_to_red', [(1, 1, 1), (1, 0, 0)], N=100)

# Define the colormaps
# cmap_white_to_coral = LinearSegmentedColormap.from_list('white_to_coral', [(1, 1, 1), (1, 0.5, 0.31)], N=100)
# cmap_white_to_teal = LinearSegmentedColormap.from_list('white_to_teal', [(1, 1, 1), (0, 0.5, 0.5)], N=100)

# Create separate colormaps for positive and negative values
norm_positive = Normalize(vmin=0.0, vmax=1)
norm_negative = Normalize(vmin=0.0, vmax=1)

sm_positive = ScalarMappable(cmap=cmap_white_to_blue, norm=norm_positive)
sm_negative = ScalarMappable(cmap=cmap_white_to_red, norm=norm_negative)

# Add colorbars underneath both ax1 and ax2
cbar_positive = fig.colorbar(sm_positive, ax=[ax1, ax2], orientation='horizontal', fraction=0.0075, pad=0.05, anchor=(0.825, -2.0))
cbar_negative = fig.colorbar(sm_negative, ax=[ax1, ax2], orientation='horizontal', fraction=0.0075, pad=0.05, anchor=(0.825, -2.0))
cbar_positive.set_label('Normalized Log Abundances', fontsize=10)


plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=600)
logger.info(f"Figure saved")