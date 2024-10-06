

### Tandem Repeat Catalog & Variation Clusters

This repo provides a general purpose genome-wide TR catalog for genotyping TR copy numbers in short read or long read sequencing data. 
It also flags TRs that reside within variation clusters, and provides locus definitions that support more accurate sequence-level analysis of these more complex regions.
Locus definitions are provided in the formats expected by many existing TR genotyping tools. 
This project is being developed as part of a collaboration between Ben Weisburd, Egor Dolzhenko, and others. 

[Release v1.0](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v1.0) is now available in draft form. 

File names that start with `repeat_catalog_v1.hg38` contain `4,863,041` TRs and are designed for repeat copy number analysis. The following formats are provided:

<br />
<table>
<tr><td><b>File suffix</b></td><td><b>Description</b></td></tr>
<tr><td>.EH.with_annotations.json.gz</td><td>For use with <a href="https://github.com/Illumina/ExpansionHunter">ExpansionHunter</a>. It includes all annotations as extra fields that are ignored by ExpansionHunter.</td></tr>
<tr><td>.EH.json.gz</td><td>For use with <a href="https://github.com/Illumina/ExpansionHunter">ExpansionHunter</a></td></tr>
<tr><td>.TRGT.bed.gz</td><td>For use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
<tr><td>.LongTR.bed.gz</td><td>For use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
<tr><td>.GangSTR.bed.gz</td><td>For use with <a href="https://github.com/gymreklab/GangSTR">GangSTR</a></td></tr>
<tr><td>.HipSTR.bed.gz</td><td>For use with <a href="https://github.com/HipSTR-Tool/HipSTR">HipSTR</a></td></tr>
<tr><td>.bed.gz</td><td>Sorted and indexed BED file for <a href="https://igv.org/">IGV</a></td></tr>
</table>
<br />

Variation clusters extend the boundaries of TRs to encompass any adjacent polymorphic regions. The exended boundaries enable more accurate sequence-level analysis, particularly in regions that contain many adjacent or multi-scale repeats with different motifs. The following files contain the extended boundary defintions:

<br />
<table>
   <tr><td><b>File name</b></td><td align="center"><b>Size</b></td><td><b>Description</b></td></tr>
   <tr><td>variation_clusters_v1.hg38.TRGT.bed.gz</td><td align="right">273,112</td><td>variation clusters for use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
   <tr><td>variation_clusters_v1.hg38.LongTR.bed.gz</td><td align="right">273,112</td><td>variation clusters for use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
   <tr><td>variation_clusters_and_isolated_repeats_v1.hg38.TRGT.bed.gz</td><td align="right">4,542,828</td><td>variation clusters + isolated repeats (ie. all tandem repeat catalog TRs that are not inside variation clusters) for use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
   <tr><td>variation_clusters_and_isolated_repeats_v1.hg38.LongTR.bed.gz</td><td align="right">4,542,828</td><td>variation clusters + isolated repeats (ie. all tandem repeat catalog TRs that are not inside variation clusters) for use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
</table>


### Goals

- Provide the catalog in the formats expected by existing TR genotyping tools for both short-read and long-read data
- Include rich annotations

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
Stats for repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz:
    4,863,041 total loci
   65,678,112 base pairs spanned by all loci (2.127% of the genome)
            0 out of  4,863,041 (  0.0%) loci define adjacent repeats
    4,863,041 total repeat intervals
    3,210,115 out of  4,863,041 ( 66.0%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
    1,567,337 out of  4,863,041 ( 32.2%) repeat intervals are homopolymers
       18,340 out of  4,863,041 (  0.4%) repeat intervals overlap each other by at least two motif lengths
           11 out of  4,863,041 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr1:82008141-82008152, chr3:78937990-78938032, chr4:1046750-1046794, chr5:52437153-52437201, chr6:34683425-34683453, chr4:107646310-107646327, chr18:52413438-52413450, chr6:150149299-150149323, chr7:40295582-40295597, chr9:35561915-35561931

Ranges:
   Motif size range: 1-833bp
   Locus size range: 1-2523bp
   Num repeats range: 1-300x repeats

   Max locus size =   2,523bp           @ chrX:71520430-71522953 (CCAGCACTTTGGGAGGCCGAGGCAGGCTGATCACTAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGGGGTGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGGGGAGGCAGGAGAATTGCCTGAACCTGAGAGGCAGAGGCTGCAGTGAGCTGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCAGTCTCAAAACAAAAAAAAAAAAAGATTTTAGTAACTTTTATCCTGTTTTAATAATACTGACTCAGAAACTATAATGTGTACTTTATAATTTACTTCCTAGATGACACTTGATTTTCTTCAAGAGCAAGATAGCTGCCCTGTGCAGTTGGTCTCCTTGAAAACTATTTTAGTTCTATCATAATTTCCTGTGATAAATATTTTGACCTTCTAAAATTTCAGAATATTGCACCAAGTAGAAAGAAAATAGGTTTTTTCTCTTTTCTTCTTCTTCCTTTTTTTTTTCTGAGAAAGAGGGAATGAGAACTTTAGTGTTCTTTCAATAGCGTTCTTATTTGTAGAAATGCATAATAGTGTCCTAGTAAGGCTTGACAATAACTCTGGTCTTCATCATATTTTGTGATAAAACTTTTGATTTAAAAAAACCTCTGATCTATTTATCATGGCAAATGGATAGAGCTTTCCTGCCTGTTTTCTTTCTTTTCTTTTTTCTTTCTTTCCTTTTTTTTCCTTTGAGCTTAGATTTTTAGAAGCACATATTTAAAAATCAGGTATAAGACTGGATGCAGTGGCTCACGCCTGTAATC)
   Min reference repeat purity   =  0.43    @ chr3:112804380-112804514 (TCT)
   Min overall mappability   =  0.00    @ chrY:56887882-56887891 (TGA)
   Base-level   purity   median: 1.000,  mean: 0.999

          chrX:    244,191 out of  4,863,041 (  5.0%) repeat intervals
          chrY:     39,257 out of  4,863,041 (  0.8%) repeat intervals
          chrM:         14 out of  4,863,041 (  0.0%) repeat intervals
   alt contigs:          0 out of  4,863,041 (  0.0%) repeat intervals

Motif size distribution:
          1bp:  1,567,337 out of  4,863,041 ( 32.2%) repeat intervals
          2bp:    978,972 out of  4,863,041 ( 20.1%) repeat intervals
          3bp:  1,432,117 out of  4,863,041 ( 29.4%) repeat intervals
          4bp:    590,787 out of  4,863,041 ( 12.1%) repeat intervals
          5bp:    177,422 out of  4,863,041 (  3.6%) repeat intervals
          6bp:     56,731 out of  4,863,041 (  1.2%) repeat intervals
       7-24bp:     43,996 out of  4,863,041 (  0.9%) repeat intervals
        25+bp:     15,679 out of  4,863,041 (  0.3%) repeat intervals

Num repeats in reference:
           1x:     10,443 out of  4,863,041 (  0.2%) repeat intervals
           2x:     38,922 out of  4,863,041 (  0.8%) repeat intervals
           3x:  1,799,189 out of  4,863,041 ( 37.0%) repeat intervals
           4x:    650,397 out of  4,863,041 ( 13.4%) repeat intervals
           5x:    356,525 out of  4,863,041 (  7.3%) repeat intervals
           6x:    151,893 out of  4,863,041 (  3.1%) repeat intervals
           7x:     85,760 out of  4,863,041 (  1.8%) repeat intervals
           8x:    257,475 out of  4,863,041 (  5.3%) repeat intervals
           9x:    352,993 out of  4,863,041 (  7.3%) repeat intervals
       10-15x:    759,188 out of  4,863,041 ( 15.6%) repeat intervals
       16-25x:    348,837 out of  4,863,041 (  7.2%) repeat intervals
       26-35x:     45,478 out of  4,863,041 (  0.9%) repeat intervals
       36-50x:      5,610 out of  4,863,041 (  0.1%) repeat intervals
         51+x:        331 out of  4,863,041 (  0.0%) repeat intervals

Reference repeat purity distribution:
          0.0:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.1:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.2:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.3:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.4:          3 out of  4,863,041 (  0.0%) repeat intervals
          0.5:         14 out of  4,863,041 (  0.0%) repeat intervals
          0.6:         44 out of  4,863,041 (  0.0%) repeat intervals
          0.7:      1,570 out of  4,863,041 (  0.0%) repeat intervals
          0.8:     12,336 out of  4,863,041 (  0.3%) repeat intervals
          0.9:     21,126 out of  4,863,041 (  0.4%) repeat intervals
          1.0:  4,827,948 out of  4,863,041 ( 99.3%) repeat intervals

Mappability distribution:
          0.0:    154,279 out of  4,863,041 (  3.2%) loci
          0.1:    214,471 out of  4,863,041 (  4.4%) loci
          0.2:    246,877 out of  4,863,041 (  5.1%) loci
          0.3:    236,856 out of  4,863,041 (  4.9%) loci
          0.4:    391,388 out of  4,863,041 (  8.0%) loci
          0.5:    561,639 out of  4,863,041 ( 11.5%) loci
          0.6:    352,273 out of  4,863,041 (  7.2%) loci
          0.7:    306,208 out of  4,863,041 (  6.3%) loci
          0.8:    337,715 out of  4,863,041 (  6.9%) loci
          0.9:    626,047 out of  4,863,041 ( 12.9%) loci
          1.0:  1,435,288 out of  4,863,041 ( 29.5%) loci

Locus sizes at each motif size:
     1bp motifs: locus size range:      1 bp to      90 bp  (median:   11 bp) based on  1,567,337 loci. Mean base purity: 1.00.  Mean mappability: 0.66
     2bp motifs: locus size range:      2 bp to     600 bp  (median:   10 bp) based on    978,972 loci. Mean base purity: 1.00.  Mean mappability: 0.76
     3bp motifs: locus size range:      3 bp to     632 bp  (median:    9 bp) based on  1,432,117 loci. Mean base purity: 1.00.  Mean mappability: 0.75
     4bp motifs: locus size range:      4 bp to     533 bp  (median:   14 bp) based on    590,787 loci. Mean base purity: 1.00.  Mean mappability: 0.68
     5bp motifs: locus size range:      5 bp to     400 bp  (median:   18 bp) based on    177,422 loci. Mean base purity: 1.00.  Mean mappability: 0.61
     6bp motifs: locus size range:      6 bp to   1,103 bp  (median:   20 bp) based on     56,731 loci. Mean base purity: 1.00.  Mean mappability: 0.62
     7bp motifs: locus size range:      7 bp to     151 bp  (median:   22 bp) based on     15,083 loci. Mean base purity: 1.00.  Mean mappability: 0.58
     8bp motifs: locus size range:      8 bp to     312 bp  (median:   25 bp) based on      7,107 loci. Mean base purity: 1.00.  Mean mappability: 0.57
     9bp motifs: locus size range:      9 bp to     153 bp  (median:   28 bp) based on      3,231 loci. Mean base purity: 1.00.  Mean mappability: 0.51
    10bp motifs: locus size range:     10 bp to     150 bp  (median:   31 bp) based on      2,713 loci. Mean base purity: 1.00.  Mean mappability: 0.50
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.log)]



