# Creating a series of bar plots, one for each row in the table
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys
import seaborn as sns

os.chdir("/Users/weisburd/p1_iter3/2023_06_23__generate_novel_STR_catalog_for_gnomAD")
args = [a for a in sys.argv if not a.startswith("--")]
catalog_stats_table_path = args[-1] if len(args) > 1 else "variant_catalog_stats.12_catalogs.plot_sort_order.tsv"
print(f"Generating plots for {catalog_stats_table_path}")

df = pd.read_table(catalog_stats_table_path)
df = df[~df["Catalog"].str.startswith("combined_catalog")]
df["Catalog"] = df["Catalog"].replace({
    "combined.51_samples.positive_loci.annotated_and_filtered.json.gz": "TR Truth Set from 51 HPRC samples",
    "hg38.hipstr_reference.adjusted.annotated_and_filtered.json.gz": "HipSTR",
    "hg38_ver17.adjusted.annotated_and_filtered.json.gz": "GangSTR v17",
    "illumina_variant_catalog.sorted.annotated_and_filtered.json.gz": "Illumina",
    "popstr_catalog_v2.annotated_and_filtered.json.gz": "popSTR",
    "repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_6bp.annotated_and_filtered.json.gz": "TRF: pure & interrupted repeats ≥ 6bp in hg38",
    "repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_6bp.annotated_and_filtered.json.gz": "TRF: only pure repeats & homopolymers ≥ 6bp in hg38",
    "repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_9bp.annotated_and_filtered.json.gz": "TRF: pure & interrupted repeats ≥ 9bp in hg38",
    "repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_9bp.annotated_and_filtered.json.gz": "TRF: only pure repeats & homopolymers ≥ 9bp in hg38,",
    "trgt_repeat_catalog.hg38.reformatted_to_motif_only.annotated_and_filtered.json.gz": "TRGT",
    "variant_catalog_without_offtargets.GRCh38.annotated_and_filtered.json.gz": "Known disease-assocaited TRs",
    "adotto_tr_catalog_v1.2.annotated_and_filtered.json.gz": "Adotto v1.2",
})

df["Catalog"] = df["Catalog"] + " (" + df["Total"] + " repeats)"

# Defining columns of interest
columns_of_interest = [
    "Homopolymers",
    "2bp motifs",
    "3bp motifs",
    "4bp motifs",
    "5bp motifs",
    "6bp motifs",
    "7+bp motifs"
]


# Number of rows in the data
n_rows = df.shape[0]

# Adjusting the plot titles and removing the "bp motifs" suffix from each x-axis label
# high resolution figure
fig, all_axes = plt.subplots(3, 4, figsize=(32, 20), sharey=False, sharex=False, constrained_layout=True, dpi=600)
df.index = range(n_rows)
i = 0
for idx_j in range(4):
    for idx_i in range(3):
        axes = all_axes[idx_i][idx_j]
        if idx_i == 0 and idx_j > 1:
            # hide the plot entirely
            axes.axis("off")
            continue

        print(df.loc[i, 'Catalog'])

        # Extracting row values and processing the data
        row_values_i = df.loc[i, columns_of_interest]
        row_values_i = row_values_i.str.rstrip("%").astype("float")
        row_values_seaborn_i = row_values_i.reset_index()
        row_values_seaborn_i.columns = ["Category", "Value"]
        row_values_seaborn_i["Category"] = [
            v.replace(" motifs", "").replace("Homopolymers", "1bp motifs") for v in row_values_seaborn_i["Category"]
        ]

        # Add a second y-axis with the absolute number of repeats

        # Plotting the bar plot for each row with updated titles and x-axis labels
        sns.barplot(x="Category", y="Value", data=row_values_seaborn_i, ax=axes, color="cornflowerblue")
        axes.set_title(f"{df.loc[i, 'Catalog']}", fontsize=14, pad=20)
        axes.set_ylabel("% of repeats", fontsize=14)
        axes.set_xlabel(" ", fontsize=14)
        axes.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")
        axes.set_ylim(0, 100)

        # add second y axis with repeat counts
        #ax2 = axes.twinx()
        #ax2.set_ylabel("Total loci", fontsize=14)
        #ax2.set_ylim(0, 100)

        # add y padding to the plot title
        axes.title.set_position([.5, 1.05])


        # add label to each bar
        for p in axes.patches:
            axes.annotate(f"{p.get_height():0.0f}%", (p.get_x() + p.get_width() / 2., p.get_height()),
                             ha="center", va="center", fontsize=12, color="black", xytext=(0, 10),
                             textcoords="offset points")
        axes.set_yticks(range(0, 101, 10))
        axes.tick_params(axis="x", labelsize=12, rotation=0, length=0)  # No x-axis tick marks
        axes.tick_params(axis="y", labelsize=12, length=0)
        i += 1
        if i >= len(df):
            break

plt.tight_layout()

# save the image
output_path = f"variant_catalog_stats.{len(df)}_catalogs.motif_sizes.png"
print(f"Saving plot to {output_path}")
plt.savefig(output_path, dpi=300)

#plt.show()
print("Done")



#%%