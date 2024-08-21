

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
<tr><td>.TRGT.bed.gz</td><td>For use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
<tr><td>.LongTR.bed.gz</td><td>For use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
<tr><td>.GangSTR.bed.gz</td><td>For use with <a href="https://github.com/gymreklab/GangSTR">GangSTR</a></td></tr>
<tr><td>.HipSTR.bed.gz</td><td>For use with <a href="https://github.com/HipSTR-Tool/HipSTR">HipSTR</a></td></tr>
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
Stats for simple_repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz:
    4,863,043 total loci
   65,678,129 base pairs spanned by all loci (2.127% of the genome)
            0 out of  4,863,043 (  0.0%) loci define adjacent repeats
    4,863,043 total repeat intervals
    3,210,117 out of  4,863,043 ( 66.0%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
    1,567,337 out of  4,863,043 ( 32.2%) repeat intervals are homopolymers
       18,340 out of  4,863,043 (  0.4%) repeat intervals overlap each other by at least two motif lengths
           11 out of  4,863,043 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr6:157039646-157039697, chr6:136309907-136309923, chr16:22663061-22663108, chr10:118053197-118053242, chr9:72024653-72024669, chr1:193880089-193880107, chr11:68644599-68644607, chr8:112342983-112342999, chr12:100569922-100569982, chr16:53171861-53171902

Ranges:
   Motif size range: 1-833bp
   Locus size range: 1-2523bp
   Num repeats range: 1-300x repeats

   Max locus size =   2,523bp           @ chrX:71520430-71522953 (CCAGCACTTTGGGAGGCCGAGGCAGGCTGATCACTAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGGGGTGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGGGGAGGCAGGAGAATTGCCTGAACCTGAGAGGCAGAGGCTGCAGTGAGCTGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCAGTCTCAAAACAAAAAAAAAAAAAGATTTTAGTAACTTTTATCCTGTTTTAATAATACTGACTCAGAAACTATAATGTGTACTTTATAATTTACTTCCTAGATGACACTTGATTTTCTTCAAGAGCAAGATAGCTGCCCTGTGCAGTTGGTCTCCTTGAAAACTATTTTAGTTCTATCATAATTTCCTGTGATAAATATTTTGACCTTCTAAAATTTCAGAATATTGCACCAAGTAGAAAGAAAATAGGTTTTTTCTCTTTTCTTCTTCTTCCTTTTTTTTTTCTGAGAAAGAGGGAATGAGAACTTTAGTGTTCTTTCAATAGCGTTCTTATTTGTAGAAATGCATAATAGTGTCCTAGTAAGGCTTGACAATAACTCTGGTCTTCATCATATTTTGTGATAAAACTTTTGATTTAAAAAAACCTCTGATCTATTTATCATGGCAAATGGATAGAGCTTTCCTGCCTGTTTTCTTTCTTTTCTTTTTTCTTTCTTTCCTTTTTTTTCCTTTGAGCTTAGATTTTTAGAAGCACATATTTAAAAATCAGGTATAAGACTGGATGCAGTGGCTCACGCCTGTAATC)
   Min fraction pure bases   =  0.43    @ chr3:112804380-112804514 (TCT)
   Min fraction pure repeats =  0.00    @ chr4:41745972-41746032 (GCN)
   Min overall mappability   =  0.00    @ chrY:56887882-56887891 (TGA)

          chrX:    244,192 out of  4,863,043 (  5.0%) repeat intervals
          chrY:     39,257 out of  4,863,043 (  0.8%) repeat intervals
          chrM:         14 out of  4,863,043 (  0.0%) repeat intervals
   alt contigs:          0 out of  4,863,043 (  0.0%) repeat intervals

Motif size distribution:
          1bp:  1,567,337 out of  4,863,043 ( 32.2%) repeat intervals
          2bp:    978,973 out of  4,863,043 ( 20.1%) repeat intervals
          3bp:  1,432,118 out of  4,863,043 ( 29.4%) repeat intervals
          4bp:    590,787 out of  4,863,043 ( 12.1%) repeat intervals
          5bp:    177,422 out of  4,863,043 (  3.6%) repeat intervals
          6bp:     56,731 out of  4,863,043 (  1.2%) repeat intervals
       7-24bp:     43,996 out of  4,863,043 (  0.9%) repeat intervals
        25+bp:     15,679 out of  4,863,043 (  0.3%) repeat intervals

Num repeats in reference:
           1x:     10,444 out of  4,863,043 (  0.2%) repeat intervals
           2x:     38,922 out of  4,863,043 (  0.8%) repeat intervals
           3x:  1,799,189 out of  4,863,043 ( 37.0%) repeat intervals
           4x:    650,397 out of  4,863,043 ( 13.4%) repeat intervals
           5x:    356,525 out of  4,863,043 (  7.3%) repeat intervals
           6x:    151,893 out of  4,863,043 (  3.1%) repeat intervals
           7x:     85,761 out of  4,863,043 (  1.8%) repeat intervals
           8x:    257,475 out of  4,863,043 (  5.3%) repeat intervals
           9x:    352,993 out of  4,863,043 (  7.3%) repeat intervals
       10-15x:    759,188 out of  4,863,043 ( 15.6%) repeat intervals
       16-25x:    348,837 out of  4,863,043 (  7.2%) repeat intervals
       26-35x:     45,478 out of  4,863,043 (  0.9%) repeat intervals
       36-50x:      5,610 out of  4,863,043 (  0.1%) repeat intervals
         51+x:        331 out of  4,863,043 (  0.0%) repeat intervals

Fraction pure bases distribution:
          0.0:        211 out of  4,863,043 (  0.0%) repeat intervals
          0.1:      1,200 out of  4,863,043 (  0.0%) repeat intervals
          0.2:      1,842 out of  4,863,043 (  0.0%) repeat intervals
          0.3:      1,711 out of  4,863,043 (  0.0%) repeat intervals
          0.4:      1,112 out of  4,863,043 (  0.0%) repeat intervals
          0.5:      5,536 out of  4,863,043 (  0.1%) repeat intervals
          0.6:     11,194 out of  4,863,043 (  0.2%) repeat intervals
          0.7:      3,902 out of  4,863,043 (  0.1%) repeat intervals
          0.8:      5,041 out of  4,863,043 (  0.1%) repeat intervals
          0.9:      3,330 out of  4,863,043 (  0.1%) repeat intervals
          1.0:  4,827,964 out of  4,863,043 ( 99.3%) repeat intervals

Mappability distribution:
          0.0:    154,279 out of  4,863,043 (  3.2%) loci
          0.1:    214,471 out of  4,863,043 (  4.4%) loci
          0.2:    246,877 out of  4,863,043 (  5.1%) loci
          0.3:    236,856 out of  4,863,043 (  4.9%) loci
          0.4:    391,389 out of  4,863,043 (  8.0%) loci
          0.5:    561,639 out of  4,863,043 ( 11.5%) loci
          0.6:    352,273 out of  4,863,043 (  7.2%) loci
          0.7:    306,208 out of  4,863,043 (  6.3%) loci
          0.8:    337,715 out of  4,863,043 (  6.9%) loci
          0.9:    626,048 out of  4,863,043 ( 12.9%) loci
          1.0:  1,435,288 out of  4,863,043 ( 29.5%) loci

Locus sizes at each motif size:
     1bp motifs: locus size range:      1 bp to      90 bp  (median:   11 bp) based on  1,567,337 loci. Mean base purity: 1.00, mean repeat purity: 1.00.  Mean mappability: 0.66
     2bp motifs: locus size range:      2 bp to     600 bp  (median:   10 bp) based on    978,973 loci. Mean base purity: 1.00, mean repeat purity: 1.00.  Mean mappability: 0.76
     3bp motifs: locus size range:      3 bp to     632 bp  (median:    9 bp) based on  1,432,118 loci. Mean base purity: 1.00, mean repeat purity: 1.00.  Mean mappability: 0.75
     4bp motifs: locus size range:      4 bp to     533 bp  (median:   14 bp) based on    590,787 loci. Mean base purity: 1.00, mean repeat purity: 0.99.  Mean mappability: 0.68
     5bp motifs: locus size range:      5 bp to     400 bp  (median:   18 bp) based on    177,422 loci. Mean base purity: 1.00, mean repeat purity: 1.00.  Mean mappability: 0.61
     6bp motifs: locus size range:      6 bp to   1,103 bp  (median:   20 bp) based on     56,731 loci. Mean base purity: 1.00, mean repeat purity: 0.98.  Mean mappability: 0.62
     7bp motifs: locus size range:      7 bp to     151 bp  (median:   22 bp) based on     15,083 loci. Mean base purity: 1.00, mean repeat purity: 0.99.  Mean mappability: 0.58
     8bp motifs: locus size range:      8 bp to     312 bp  (median:   25 bp) based on      7,107 loci. Mean base purity: 1.00, mean repeat purity: 0.99.  Mean mappability: 0.57
     9bp motifs: locus size range:      9 bp to     153 bp  (median:   28 bp) based on      3,231 loci. Mean base purity: 1.00, mean repeat purity: 0.97.  Mean mappability: 0.51
    10bp motifs: locus size range:     10 bp to     150 bp  (median:   31 bp) based on      2,713 loci. Mean base purity: 1.00, mean repeat purity: 0.96.  Mean mappability: 0.50
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.merge_and_annotate_loci.log)]



