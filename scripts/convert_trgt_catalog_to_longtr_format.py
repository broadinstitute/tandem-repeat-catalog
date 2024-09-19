"""This script takes a TRGT catalog BED file and outputs a new BED file in LongTR format."""

import argparse
import collections
import gzip
import ijson
import os
import simplejson as json
import re
import sys
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog import convert_expansion_hunter_record_to_trgt_row

def compute_dominant_repeat_unit_length(info_fields_dict, known_pathogenic_reference_regions_lookup):
    """Compute the dominant repeat unit length for a TR locus.

    Args:
        info_fields_dict (dict): A dictionary of info fields for a TR locus
        known_pathogenic_reference_regions_lookup (dict): A dictionary mapping locus IDs to reference regions for all
            known disease-associated loci

    Return:
        int: The dominant repeat unit length for the TR locus.
    """
    repeat_unit_lengths = []
    for locus_id in info_fields_dict["ID"].split(","):
        if locus_id in known_pathogenic_reference_regions_lookup:
            # get the repeat unit length
            _, repeat_unit_length = max((len(ru), len(ru)) for ru in info_fields_dict["MOTIFS"].split(","))
            repeat_unit_lengths.append((10**9, repeat_unit_length))
        elif locus_id.count("-") == 3:

            chrom, start_0based, end, repeat_unit = locus_id.split("-")
            repeat_unit_lengths.append((int(end) - int(start_0based), len(repeat_unit)))
        else:
            raise ValueError(f"Unexpected locus_id '{locus_id}'")

    _, repeat_unit_length = max(repeat_unit_lengths)
    return repeat_unit_length


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
                        "containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
                        "these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
    parser.add_argument("-o", "--output-bed-path",
                        help="Path of output BED file.")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("input_trgt_catalog_bed_path", help="Path of the input TRGT catalog BED file")
    args = parser.parse_args()

    if ".bed" not in args.input_trgt_catalog_bed_path:
        parser.error(f"Input path {args.input_trgt_catalog_bed_path} path must does not have a '.bed' suffix")

    if not args.output_bed_path:
        args.output_bed_path = re.sub(".bed(.b?gz)?$", "", args.input_trgt_catalog_bed_path)
        args.output_bed_path = args.output_bed_path.replace(".TRGT", "")
        args.output_bed_path += ".LongTR.bed"
    elif not args.output_bed_path.endswith(".bed"):
        parser.error("--output-bed-path must have a '.bed' suffix")

    fopen = gzip.open if args.known_pathogenic_loci_json_path.endswith("gz") else open
    with fopen(args.known_pathogenic_loci_json_path, "rt") as f:
        known_pathogenic_loci = json.load(f)
        known_pathogenic_reference_regions_lookup = {}
        for locus in known_pathogenic_loci:
            if isinstance(locus["ReferenceRegion"], list):
                assert isinstance(locus["VariantId"], list)
                assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
                for variant_id, reference_region in zip(locus["VariantId"], locus["ReferenceRegion"]):
                    known_pathogenic_reference_regions_lookup[variant_id] = reference_region
            else:
                known_pathogenic_reference_regions_lookup[locus["LocusId"]] = locus["ReferenceRegion"]

    # get locus ids of TRs in variation clusters
    counter = collections.Counter()
    output_bed_file = open(args.output_bed_path, "wt")
    fopen = gzip.open if args.input_trgt_catalog_bed_path.endswith("gz") else open
    with fopen(args.input_trgt_catalog_bed_path, "rt") as f:
        if args.show_progress_bar:
            f = tqdm.tqdm(f, unit=" records", unit_scale=True)

        for line in f:
            counter["total"] += 1
            fields = line.strip("\n").split("\t")
            chrom = fields[0]
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            info_fields = fields[3]

            info_fields_dict = {}
            for key_value in info_fields.split(";"):
                key_value = key_value.split("=")
                if len(key_value) != 2:
                    print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
                    continue
                key, value = key_value
                info_fields_dict[key] = value

            repeat_unit_length = compute_dominant_repeat_unit_length(
                info_fields_dict, known_pathogenic_reference_regions_lookup)

            output_bed_file.write("\t".join(map(str, [
                chrom,
                start_0based + 1,  # LongTR BED files use 1-based coords.
                end_1based,
                repeat_unit_length,
                int((end_1based - start_0based)/repeat_unit_length),
                info_fields_dict["ID"],
            ])) + "\n")

    os.system(f"bgzip -f {args.output_bed_path}")

    print(f"Wrote {counter['total']:,d} rows to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()

