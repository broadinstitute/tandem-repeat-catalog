

### Simple Repeat Catalog

This repo provides a general purpose genome-wide TR catalog for genotyping TRs in short read or long read sequencing data. 
It is being developed as part of a collaboration between Ben Weisburd, Egor Dolzhenko, and others. 

[Release v1.0](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v1.0) is now available in draft form. 

It includes filenames that start with `simple_repeat_catalog_v1.hg38` and have the following suffixes:

<br />
<table>
<tr><td><b>File suffix</b></td><td><b>Description</b></td></tr>
<tr><td>.merged_and_annotated.json.gz</td><td>For use with <a href="https://github.com/Illumina/ExpansionHunter">ExpansionHunter</a>. It includes all annotations as extra fields that are ignored by ExpansionHunter.</td></tr>
<tr><td>.bed.gz</td><td>Sorted and indexed BED file for IGV</td></tr>
<tr><td>.TRGT.bed</td><td>For use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
<tr><td>.LongTR.bed</td><td>For use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
<tr><td>.GangSTR.bed</td><td>For use with <a href="https://github.com/gymreklab/GangSTR">GangSTR</a></td></tr>
<tr><td>.HipSTR.bed</td><td>For use with <a href="https://github.com/HipSTR-Tool/HipSTR">HipSTR</a></td></tr>
</table>
<br />

Additionally, `variation_clusters_v1.hg38.TRGT.bed.gz` contains variation clusters in a format that can be passed to TRGT to genotype wider regions around the simple repeats defined in the catalog above. 

### Goals

- Create a catalog that is as sensitive (includes all polymorphic TR loci in the human genome) and specific (excludes non-polymorphic loci) as possible
- Share the catalog in formats that can be used directly with most TR genotyping tools for both short read and long read data. 
- Provide rich annotations

Although our initial focus is on the human genome, we'd also love to extend this work to plants and other species. Please consider creating a GitHub issue or reaching out by email if you are interested in this.

### Background

Tandem repeats (TRs) are regions of the genome that consist of consecutive copies of some motif sequence. For example, `CAGCAGCAG` is a tandem repeat of the `CAG` motif. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs, which specify the genomic start and end coordinates of each tandem repeat region, as well as the one or more motifs that repeat there. 

For example, if a hypothetical region at the beginning of `chrX` consisted of the following nucleotide sequence:  
`ATCAGTAGA ATATATATAT CAGACAGCAGCAG TGAGTGCGTAC...`  
it could be represented in a repeat catalog as two entries:  
`chrX:10-19 (AT)*`  
`chrX:20-32 (CAG)*`   
indicating that a repeat of the `AT` motif occurs between positions 10 and 19 (inclusive), and of the `CAG` motif between positions 20 and 32.
A genome-wide catalog would contain such entries for all repeat regions of interest found anywhere in the genome. 


### Catalog Stats

The following catalog stats for v1.0 were computed using [str_analysis/compute_catalog_stats.py](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/compute_catalog_stats.py):

```
Stats for simple_repeat_catalog_v1.hg38.1_to_1000bp_motifs.merged_and_annotated.json.gz:
    4,623,360 total loci
   64,363,319 base pairs spanned by all loci (2.084% of the genome)
            0 out of  4,623,360 (  0.0%) loci define adjacent repeats
    4,623,360 total repeat intervals
    3,005,968 out of  4,623,360 ( 65.0%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
    1,395,502 out of  4,623,360 ( 30.2%) repeat intervals are homopolymers
       17,637 out of  4,623,360 (  0.4%) repeat intervals overlap each other by at least two motif lengths
           11 out of  4,623,360 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr4:10997247-10997259, chr17:54429149-54429165, chr7:75324911-75324935, chr10:34074010-34074030, chr1:17191894-17191906, chr6:1512741-1512765, chr14:28801139-28801165, chr6:12921458-12921494, chr6:144337148-144337191, ch
r1:207100053-207100097

Ranges:
   Motif size range: 1-833bp
   Locus size range: 1-2523bp
   Num repeats range: 1-300x repeats
   Maximum locus size = 2,523bp               @ chrX:71520430-71522953

          chrX:    234,392 out of  4,623,360 (  5.1%) repeat intervals
          chrY:     38,171 out of  4,623,360 (  0.8%) repeat intervals
          chrM:         14 out of  4,623,360 (  0.0%) repeat intervals
   alt contigs:          0 out of  4,623,360 (  0.0%) repeat intervals

Motif size distribution:
          1bp:  1,395,502 out of  4,623,360 ( 30.2%) repeat intervals
          2bp:    939,879 out of  4,623,360 ( 20.3%) repeat intervals
          3bp:  1,421,178 out of  4,623,360 ( 30.7%) repeat intervals
          4bp:    578,904 out of  4,623,360 ( 12.5%) repeat intervals
          5bp:    175,107 out of  4,623,360 (  3.8%) repeat intervals
          6bp:     55,842 out of  4,623,360 (  1.2%) repeat intervals
       7-24bp:     42,231 out of  4,623,360 (  0.9%) repeat intervals
        25+bp:     14,717 out of  4,623,360 (  0.3%) repeat intervals

Num repeats in reference:
           1x:      7,236 out of  4,623,360 (  0.2%) repeat intervals
           2x:     33,399 out of  4,623,360 (  0.7%) repeat intervals
           3x:  1,775,953 out of  4,623,360 ( 38.4%) repeat intervals
           4x:    621,237 out of  4,623,360 ( 13.4%) repeat intervals
           5x:    332,203 out of  4,623,360 (  7.2%) repeat intervals
           6x:    140,297 out of  4,623,360 (  3.0%) repeat intervals
           7x:     78,379 out of  4,623,360 (  1.7%) repeat intervals
           8x:    108,624 out of  4,623,360 (  2.3%) repeat intervals
           9x:    346,452 out of  4,623,360 (  7.5%) repeat intervals
       10-15x:    756,349 out of  4,623,360 ( 16.4%) repeat intervals
       16-25x:    365,283 out of  4,623,360 (  7.9%) repeat intervals
       26-35x:     51,261 out of  4,623,360 (  1.1%) repeat intervals
       36-50x:      6,312 out of  4,623,360 (  0.1%) repeat intervals
         51+x:        375 out of  4,623,360 (  0.0%) repeat intervals

Fraction pure bases distribution:
          0.0:        176 out of  4,623,360 (  0.0%) repeat intervals
          0.1:      1,020 out of  4,623,360 (  0.0%) repeat intervals
          0.2:      1,575 out of  4,623,360 (  0.0%) repeat intervals
          0.3:      1,515 out of  4,623,360 (  0.0%) repeat intervals
          0.4:        977 out of  4,623,360 (  0.0%) repeat intervals
          0.5:      4,769 out of  4,623,360 (  0.1%) repeat intervals
          0.6:      9,805 out of  4,623,360 (  0.2%) repeat intervals
          0.7:      3,718 out of  4,623,360 (  0.1%) repeat intervals
          0.8:     23,714 out of  4,623,360 (  0.5%) repeat intervals
          0.9:     54,214 out of  4,623,360 (  1.2%) repeat intervals
          1.0:  4,521,877 out of  4,623,360 ( 97.8%) repeat intervals

Mappability distribution:
          0.0:    134,451 out of  4,623,360 (  2.9%) loci
          0.1:    154,470 out of  4,623,360 (  3.3%) loci
          0.2:    189,922 out of  4,623,360 (  4.1%) loci
          0.3:    266,520 out of  4,623,360 (  5.8%) loci
          0.4:    399,582 out of  4,623,360 (  8.6%) loci
          0.5:    571,542 out of  4,623,360 ( 12.4%) loci
          0.6:    343,106 out of  4,623,360 (  7.4%) loci
          0.7:    308,318 out of  4,623,360 (  6.7%) loci
          0.8:    334,406 out of  4,623,360 (  7.2%) loci
          0.9:    646,363 out of  4,623,360 ( 14.0%) loci
          1.0:  1,274,680 out of  4,623,360 ( 27.6%) loci
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.merge_and_annotate_loci.log)]



