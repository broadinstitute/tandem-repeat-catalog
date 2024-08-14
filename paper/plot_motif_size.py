# Creating a series of bar plots, one for each row in the table
import matplotlib.pyplot as plt
import os
import pandas as pd

import seaborn as sns

catalog_stats_table_path = "compare_catalogs/combined_catalog_stats.all_13_catalogs.tsv"

print(f"Generating plots for {catalog_stats_table_path}")
df = pd.read_table(catalog_stats_table_path)
df = df[df["catalog"] != "known_disease_associated_loci"]
df["catalog"] = df["catalog"].replace({
        "primary_known_disease_associated_loci": "Known disease-associated loci",
        "mukamel_VNTR_catalog": " VNTR catalog from Mukamel et al. 2021",
        "comprehensive_catalog_from_Chiu_et_al": "Comprehensive catalog from Chiu et al. 2024",
        "perfect_repeats_in_hg38": "All perfect repeats (≥ 3x and ≥ 9bp) in hg38",
        "illumina_catalog": "Illumina catalog of polymorphic loci in 1kGP",
        "polymorphic_loci_in_HPRC_assemblies": "Polymorphic loci in 51 HPRC assemblies",
    })
df["catalog"] = df["catalog"] + " (" + df["total"].apply(lambda t: f"{t:,d}") + " loci)"

#'count_1bp_motifs', 'count_2bp_motifs', 'count_3bp motifs', 'count_4bp_motifs',
#      'count_5bp_motifs', 'count_6bp_motifs', 'count_7+bp_motifs',

df.rename(columns={
    "count_1bp_motifs": "1bp motifs",
    "count_2bp_motifs": "2bp motifs",
    "count_3bp motifs": "3bp motifs",
    "count_4bp_motifs": "4bp motifs",
    "count_5bp_motifs": "5bp motifs",
    "count_6bp_motifs": "6bp motifs",
    "count_7+bp_motifs": "7+bp motifs",
    #"count_7-24bp_motifs": "7-24bp motifs",
    #"count_25+bp_motifs": "25+bp motifs",
}, inplace=True)


#print(df.columns)

# Defining columns of interest
columns_of_interest = [
    "1bp motifs",
    "2bp motifs",
    "3bp motifs",
    "4bp motifs",
    "5bp motifs",
    "6bp motifs",
    "7+bp motifs",
    #"7-24bp motifs",
    #"25+bp motifs",
]

df.index = range(len(df))

fig, all_axes = plt.subplots(4, 3, figsize=(20, 28), sharey=True, sharex=False, constrained_layout=True, dpi=600)
i = 0
for idx_i in range(4):
    for idx_j in range(3):
        axes = all_axes[idx_i][idx_j]

        current_catalog_name = df.loc[i, 'catalog'].replace("_", " ")
        total = int(df.loc[i, "total"])

        print(current_catalog_name)

        row_values_i = df.loc[i, columns_of_interest]
        row_values_seaborn_i = row_values_i.reset_index()
        row_values_seaborn_i.columns = ["Category", "Count"]
        row_values_seaborn_i["Category"] = [
            v.replace(" motifs", "") for v in row_values_seaborn_i["Category"]
        ]
        row_values_seaborn_i["Value"] = 100 * row_values_seaborn_i["Count"] / total
        print(row_values_seaborn_i)
        # Add a second y-axis with the absolute number of repeats

        # Plotting the bar plot for each row with updated titles and x-axis labels
        sns.barplot(x="Category", y="Value", data=row_values_seaborn_i, ax=axes, color="cornflowerblue")
        axes.set_title(f"{current_catalog_name}", fontsize=15, pad=20)
        axes.set_ylabel("% of loci", fontsize=15)
        axes.set_xlabel(" ", fontsize=15)
        axes.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")
        axes.set_ylim(0, 100)

        # add second y axis with repeat counts
        #ax2 = axes.twinx()
        #ax2.set_ylabel("Total loci", fontsize=14)
        #ax2.set_ylim(0, 100)

        # add y padding to the plot title
        axes.title.set_position([.5, 1.05])


        # add label to each bar
        for locus_count, p in zip(row_values_seaborn_i["Count"], axes.patches):
            axes.annotate(f"{int(locus_count):,d}", (p.get_x() + p.get_width() / 2., p.get_height()),
                             ha="center", va="center", fontsize=11, color="black", xytext=(0, 10),
                             textcoords="offset points")

        axes.set_yticks(range(0, 101, 10))
        axes.tick_params(axis="x", labelsize=14, rotation=0, length=0)  # No x-axis tick marks
        axes.tick_params(axis="y", labelsize=14, length=0)
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