import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
from pprint import pprint
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--stats-table-path", default="compare_catalogs/combined_catalog_stats.all_13_catalogs.tsv")

    parser.add_argument("--skip-motif-size-plots", action="store_true", help="Skip plotting motif size histograms")
    parser.add_argument("--skip-locus-size-plots", action="store_true", help="Skip plotting locus size histograms")
    parser.add_argument("--only-plot-motif-size", action="store_true", help="Only plot motif size histograms")
    parser.add_argument("--only-plot-locus-size", action="store_true", help="Only plot locus size histograms")

    parser.add_argument("--grid-width", type=int, default=3, help="Number of columns in the trellis grid")
    parser.add_argument("--grid-height", type=int, default=4, help="Number of rows in the trellis grid")

    args = parser.parse_args()

    catalog_stats_table_path = args.stats_table_path
    if not os.path.isfile(catalog_stats_table_path):
        parser.error(f"File not found: {catalog_stats_table_path}")

    if args.grid_width * args.grid_height != 12:
        parser.error("grid_width * grid_height must equal 12")

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

    print(f"Parsed {len(df)} catalogs from {catalog_stats_table_path}")

    #pprint(list(df.columns))
    output_paths = []
    for histogram_type in "motif_sizes", "locus_sizes":
        output_path = create_plot(histogram_type, df, args)
        output_paths.append(output_path)

    print("Done generating", " and ".join(output_paths))


def create_plot(histogram_type, df, args):
    grid_rows = args.grid_height
    grid_columns = args.grid_width

    output_path = f"variant_catalog_stats.{len(df)}_catalogs.{grid_rows}x{grid_columns}.{histogram_type}.png"
    print(f"Generating plots for {output_path}")
    columns_of_interest = []
    column_name_map = {}
    if histogram_type == "motif_sizes":
        for label in list(range(1, 7)) + ["7-24", "25+"]:
            new_column_name = f"{label}bp motifs"
            column_name_map[f"count_{label}bp_motifs"] = new_column_name
            columns_of_interest.append(new_column_name)
    elif histogram_type == "locus_sizes":
        for label in list(range(1, 25)) + ["25-50", "51+"]:
            new_column_name = f"{label}x"
            column_name_map[f"motif_sizes_per_locus_size:{label}x"] = new_column_name
            columns_of_interest.append(new_column_name)
    else:
        raise ValueError(f"Invalid histogram_type: {histogram_type}")
    df.rename(columns=column_name_map, inplace=True)
    df.index = range(len(df))

    fig, all_axes = plt.subplots(grid_rows, grid_columns, figsize=(7 * grid_columns, 7 * grid_rows), sharey=True,
                                 sharex=False, constrained_layout=True)  # , dpi=600)
    i = 0
    for idx_i in range(grid_rows):
        for idx_j in range(grid_columns):
            axes = all_axes[idx_i][idx_j]

            # if idx_i >= 1 or idx_j >= 1:
            #    # don't plot axes
            #    axes.axis("off")
            #    continue

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

            sns.barplot(x="Category", y="Value", data=row_values_seaborn_i, ax=axes, color="cornflowerblue")

            axes.set_title(f"{current_catalog_name}\n\n{total:,d} loci", fontsize=20, pad=20)
            axes.set_ylabel("% of loci", fontsize=18)
            axes.set_xlabel(" ", fontsize=18)
            axes.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")

            # add second y axis with repeat counts
            # ax2 = axes.twinx()
            # ax2.set_ylabel("Total loci", fontsize=14)
            # ax2.set_ylim(0, 100)

            # remote spines
            axes.spines["top"].set_visible(False)
            axes.spines["right"].set_visible(False)
            #if idx_j > 0:
            #    axes.spines["left"].set_visible(False)

            # add y padding to the plot title
            axes.title.set_position([.5, 1.05])

            axes.tick_params(axis="x", labelsize=14, rotation=0, length=0)  # No x-axis tick marks
            axes.tick_params(axis="y", labelsize=14, length=0)

            if histogram_type == "motif_sizes":
                y_max = 100
                counts_and_patches = zip(row_values_seaborn_i["Count"], axes.patches)
            elif histogram_type == "locus_sizes":
                y_max = 70
                largest_5_values = set(list(sorted(row_values_seaborn_i["Count"]))[-5:])
                counts_and_patches = [
                    (locus_count, p) for locus_count, p in zip(row_values_seaborn_i["Count"], axes.patches)
                    if locus_count == 0 or locus_count in largest_5_values
                ]
            else:
                raise ValueError(f"Invalid histogram_type: {histogram_type}")

            axes.set_ylim(0, y_max)
            axes.set_yticks(range(0, y_max + 1, 10))

            for locus_count, p in counts_and_patches:
                x_pos = p.get_x() + 0.5
                y_pos = p.get_height()
                axes.annotate(f"{int(locus_count):,d}", (x_pos, y_pos),
                              ha="center", va="bottom", fontsize=15, color="gray", xytext=(0, 10),
                              textcoords="offset points", rotation=90)

            axes.set_xticklabels(axes.get_xticklabels(), rotation=90)
            axes.set_xlabel(axes.get_xlabel(), fontsize=18)

            i += 1
            if i >= len(df):
                break

    plt.tight_layout()

    print(f"Saving plot to {output_path}")
    plt.savefig(output_path)  # , dpi=300)

    return output_path


if __name__ == "__main__":
    main()