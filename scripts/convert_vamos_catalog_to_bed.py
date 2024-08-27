import os
import pandas as pd
from pprint import pprint

if os.getcwd().endswith("str-truth-set-v2"):
	os.chdir("str-truth-set/ref/other")
elif os.getcwd().endswith("str-truth-set"):
	os.chdir("ref/other")

# source: https://github.com/ChaissonLab/vamos
version = "v2.1"
df = pd.read_table(
	f"https://zenodo.org/records/11625069/files/vamos.motif.hg38.{version}.e0.1.tsv.gz?download=1",
	compression="gzip",
	names=["chrom", "start_1based", "end_1based", "motifs", "version", "STR_or_VNTR", "motif_size", "score1", "score2", "segdup", "is_coding"],
)
df["start_0based"] = df["start_1based"] - 1
df["locus_size"] = df["end_1based"] - df["start_0based"]
df["motif1"] = df["motifs"].str.split(",").str[0]
df["motif_count"] = df["motifs"].str.count(",") + 1
df["is_coding"] = df["is_coding"] == "coding"
df["segdup"] = df["segdup"] == "segDup"

output_path_prefix = f"vamos_catalog.{version}"
df.to_csv(f"{output_path_prefix}.tsv.gz", sep="\t", index=False, header=True)

df = df[["chrom", "start_0based", "end_1based", "motif1", "motif_size"]]
df.to_csv(f"{output_path_prefix}.bed", sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_path_prefix}.bed")
os.system(f"tabix -f {output_path_prefix}.bed.gz")
print(f"Wrote {len(df):,d} loci to {output_path_prefix}.bed.gz")

notes = """Notes: 
len(df) == 1,175,953
df["version"] is always "v-2.0"
df["STR_or_VNTR"] is:
	STR     805485
	VNTR    370468
	
motif sizes:
	1      175625
	2      194556
	3	  	76515
	4      175497
	5       86683
	6       96609
	7       56180
	8        9294
	9	     5361
        ...
  491         1  
  492		  2
  494         1
  496         1

"""

pprint(df.iloc[0])
