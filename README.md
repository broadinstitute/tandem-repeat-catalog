

### Colab Repeat Catalog

This repo provides a general purpose genome-wide TR catalog for genotyping TRs in short read or long read sequencing data. 
It is being developed as part of a collaboration between Ben Weisburd, Egor Dolzhenko, and others. 

[Release v1.0](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v1.0) is now available in draft form. 

It includes filenames that start with `repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38` and have the following suffixes:

<br />
<table>
<tr><td><b>File suffix</b></td><td><b>Description</b></td></tr>
<tr><td>.merged_and_annotated.json.gz</td><td>For use with <a href="https://github.com/Illumina/ExpansionHunter">ExpansionHunter</a>. It includes all annotations as extra fields which ExpansionHunter ignores.</td></tr>
<tr><td>.bed.gz</td><td>Sorted and indexed BED file for IGV</td></tr>
<tr><td>.TRGT.bed</td><td>For use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
<tr><td>.LongTR.bed</td><td>For use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
<tr><td>.GangSTR.bed</td><td>For use with <a href="https://github.com/gymreklab/GangSTR">GangSTR</a></td></tr>
<tr><td>.HipSTR.bed</td><td>For use with <a href="https://github.com/HipSTR-Tool/HipSTR">HipSTR</a></td></tr>
</table>
<br />

Additionally, `variation_clusters.hg38.TRGT.bed` contains variation clusters in a format that can be passed to TRGT to genotype wider regions around the simple repeats defined in the catalog above. 

### Goals

- Create a catalog that is as sensitive (includes all polymorphic TR loci in the human genome) and specific (excludes non-polymorphic loci) as possible
- Provide rich annotations
- Share the catalog in formats that can be used directly with most TR genotyping tools for both short read and long read data. 

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
Stats for repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.merged_and_annotated.json.gz:
    3,227,828 total loci
   45,754,784 base pairs spanned by all loci (1.482% of the genome)
            0 out of  3,227,828 (  0.0%) loci define adjacent repeats
    3,227,828 total repeat intervals
    1,610,453 out of  3,227,828 ( 49.9%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
            0 out of  3,227,828 (  0.0%) repeat intervals are homopolymers
       17,061 out of  3,227,828 (  0.5%) repeat intervals overlap each other by at least two motif lengths
           11 out of  3,227,828 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr12:47354714-47354726, chr11:1543057-1543066, chr15:58276212-58276239, chr5:67345051-67345071, chr5:86960885-86960921, chr7:145020235-145020267, chrX:87717009-87717048, chr6:1827076-1827124, chr2:66602257-66602281, chr17:35469321-35469339

Ranges:
   Motif size range: 2-833bp
   Locus size range: 2-2523bp
   Num repeats range: 1-300x repeats

   Maximum locus size = 2523bp               @ chrX:71520430-71522953 (CCAGCACTTTGGGAGGCCGAGGCAGGCTGATCACTAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGGGGTGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGGGGAGGCAGGAGAATTGCCTGAACCTGAGAGGCAGAGGCTGCAGTGAGCTGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCAGTCTCAAAACAAAAAAAAAAAAAGATTTTAGTAACTTTTATCCTGTTTTAATAATACTGACTCAGAAACTATAATGTGTACTTTATAATTTACTTCCTAGATGACACTTGATTTTCTTCAAGAGCAAGATAGCTGCCCTGTGCAGTTGGTCTCCTTGAAAACTATTTTAGTTCTATCATAATTTCCTGTGATAAATATTTTGACCTTCTAAAATTTCAGAATATTGCACCAAGTAGAAAGAAAATAGGTTTTTTCTCTTTTCTTCTTCTTCCTTTTTTTTTTCTGAGAAAGAGGGAATGAGAACTTTAGTGTTCTTTCAATAGCGTTCTTATTTGTAGAAATGCATAATAGTGTCCTAGTAAGGCTTGACAATAACTCTGGTCTTCATCATATTTTGTGATAAAACTTTTGATTTAAAAAAACCTCTGATCTATTTATCATGGCAAATGGATAGAGCTTTCCTGCCTGTTTTCTTTCTTTTCTTTTTTCTTTCTTTCCTTTTTTTTCCTTTGAGCTTAGATTTTTAGAAGCACATATTTAAAAATCAGGTATAAGACTGGATGCAGTGGCTCACGCCTGTAATC)
   Minimum fraction pure bases = 0.17      @ chr8:49294435-49294500 (AGAT)
   Minimum fraction pure repeats = 0.00    @ chr4:41745972-41746032 (GCN)
   Minimum overall mappability = 0.00       @ chrY:56887882-56887891 (TGA)

          chrX:    167,667 out of  3,227,828 (  5.2%) repeat intervals
          chrY:     29,107 out of  3,227,828 (  0.9%) repeat intervals
          chrM:         14 out of  3,227,828 (  0.0%) repeat intervals
   alt contigs:          0 out of  3,227,828 (  0.0%) repeat intervals

Motif size distribution:
          1bp:          0 out of  3,227,828 (  0.0%) repeat intervals
          2bp:    939,869 out of  3,227,828 ( 29.1%) repeat intervals
          3bp:  1,421,176 out of  3,227,828 ( 44.0%) repeat intervals
          4bp:    578,898 out of  3,227,828 ( 17.9%) repeat intervals
          5bp:    175,096 out of  3,227,828 (  5.4%) repeat intervals
          6bp:     55,842 out of  3,227,828 (  1.7%) repeat intervals
       7-24bp:     42,231 out of  3,227,828 (  1.3%) repeat intervals
        25+bp:     14,716 out of  3,227,828 (  0.5%) repeat intervals

Num repeats in reference:
           1x:      6,840 out of  3,227,828 (  0.2%) repeat intervals
           2x:     33,264 out of  3,227,828 (  1.0%) repeat intervals
           3x:  1,775,635 out of  3,227,828 ( 55.0%) repeat intervals
           4x:    620,843 out of  3,227,828 ( 19.2%) repeat intervals
           5x:    331,430 out of  3,227,828 ( 10.3%) repeat intervals
           6x:    139,178 out of  3,227,828 (  4.3%) repeat intervals
           7x:     74,654 out of  3,227,828 (  2.3%) repeat intervals
           8x:     43,180 out of  3,227,828 (  1.3%) repeat intervals
           9x:     31,615 out of  3,227,828 (  1.0%) repeat intervals
       10-15x:    108,397 out of  3,227,828 (  3.4%) repeat intervals
       16-25x:     57,478 out of  3,227,828 (  1.8%) repeat intervals
       26-35x:      4,839 out of  3,227,828 (  0.1%) repeat intervals
       36-50x:        362 out of  3,227,828 (  0.0%) repeat intervals
         51+x:        113 out of  3,227,828 (  0.0%) repeat intervals

Fraction pure bases distribution:
          0.0:        176 out of  3,227,828 (  0.0%) repeat intervals
          0.1:      1,020 out of  3,227,828 (  0.0%) repeat intervals
          0.2:      1,576 out of  3,227,828 (  0.0%) repeat intervals
          0.3:      1,515 out of  3,227,828 (  0.0%) repeat intervals
          0.4:        977 out of  3,227,828 (  0.0%) repeat intervals
          0.5:      4,750 out of  3,227,828 (  0.1%) repeat intervals
          0.6:      9,768 out of  3,227,828 (  0.3%) repeat intervals
          0.7:      3,942 out of  3,227,828 (  0.1%) repeat intervals
          0.8:     24,042 out of  3,227,828 (  0.7%) repeat intervals
          0.9:     36,332 out of  3,227,828 (  1.1%) repeat intervals
          1.0:  3,143,730 out of  3,227,828 ( 97.4%) repeat intervals

Mappability distribution:
          0.0:    102,868 out of  3,227,828 (  3.2%) loci
          0.1:     99,509 out of  3,227,828 (  3.1%) loci
          0.2:    121,832 out of  3,227,828 (  3.8%) loci
          0.3:    175,517 out of  3,227,828 (  5.4%) loci
          0.4:    230,667 out of  3,227,828 (  7.1%) loci
          0.5:    309,570 out of  3,227,828 (  9.6%) loci
          0.6:    227,112 out of  3,227,828 (  7.0%) loci
          0.7:    222,580 out of  3,227,828 (  6.9%) loci
          0.8:    258,585 out of  3,227,828 (  8.0%) loci
          0.9:    513,917 out of  3,227,828 ( 15.9%) loci
          1.0:    965,671 out of  3,227,828 ( 29.9%) loci
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.merge_and_annotate_loci.log)], including:

```
Merged 58,226 out of 3,286,072 (  1.8%) adjacent loci that were within 1bp of each other and had the same motif after cyclic shift

Sources of loci in the final catalog:
          82 out of 3,286,072 ( 0.0%) Known disease-associated loci (variant_catalog_without_offtargets.GRCh38.split.filtered.json)
     173,879 out of 3,286,072 ( 5.3%) Illumina catalog of 174k polymorphic loci (illumina_variant_catalog.sorted.filtered.json)
   3,053,992 out of 3,286,072 (92.9%) All pure repeats in hg38 (hg38_repeats.motifs_2_to_1000bp.repeats_3x_and_spans_9bp.filtered.json)
      58,119 out of 3,286,072 ( 1.8%) TR variants from 51 HPRC assemblies (combined.51_samples.positive_loci.filtered.json)
```


