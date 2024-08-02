

### Colab Repeat Catalog

This repo provides a general purpose genome-wide TR catalog for genotyping TRs in short read or long read sequencing data. 

It's named *Colab Repeat Catalog* as it is a collaboration between multiple people including Ben Weisburd, Egor Dolzhenko, and others. 

[Release v1.0](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v1.0) is now in draft form. It includes filenames that start with `repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38` and have the following suffixes:

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

### Goals

- Create a catalog that is as sensitive (includes all polymorphic TR loci in the human genome) and specific (excludes non-polymorphic loci) as possible
- Provide rich annotations
- Share the catalog in formats that can be used directly with widely-used TR genotyping tools

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

The following catalog stats for v0.9 were computed using [str_analysis/compute_catalog_stats.py](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/compute_catalog_stats.py):

```
Stats for repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.merged_and_annotated.json.gz:
    3,226,630 total loci
   45,897,641 base pairs spanned by all loci (1.486% of the genome)
            0 out of  3,226,630 (  0.0%) loci define adjacent repeats
    3,226,630 total repeat intervals
    1,509,050 out of  3,226,630 ( 46.8%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
            0 out of  3,226,630 (  0.0%) repeat intervals are homopolymers
       16,341 out of  3,226,630 (  0.5%) repeat intervals overlap each other by at least two motif lengths
           11 out of  3,226,630 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr1:192003140-192003240, chr14:94601030-94601072, chr2:7489031-7489043, chr4:49113861-49113896, chr8:97735331-97735348, chr12:52918096-52918120, chr2:228457415-228457451, chr16:50476869-50476885, chr4:74413624-74413689, chr5:4438330-4438352

Ranges:
   Motif size range: 2-833bp
   Locus size range: 2-2523bp
   Num repeats range: 1-300x repeats

   Maximum locus size = 2523bp               @ chrX:71520430-71522953 (CCAGCACTTTGGGAGGCCGAGGCAGGCTGATCACTAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGGGGTGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGGGGAGGCAGGAGAATTGCCTGAACCTGAGAGGCAGAGGCTGCAGTGAGCTGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCAGTCTCAAAACAAAAAAAAAAAAAGATTTTAGTAACTTTTATCCTGTTTTAATAATACTGACTCAGAAACTATAATGTGTACTTTATAATTTACTTCCTAGATGACACTTGATTTTCTTCAAGAGCAAGATAGCTGCCCTGTGCAGTTGGTCTCCTTGAAAACTATTTTAGTTCTATCATAATTTCCTGTGATAAATATTTTGACCTTCTAAAATTTCAGAATATTGCACCAAGTAGAAAGAAAATAGGTTTTTTCTCTTTTCTTCTTCTTCCTTTTTTTTTTCTGAGAAAGAGGGAATGAGAACTTTAGTGTTCTTTCAATAGCGTTCTTATTTGTAGAAATGCATAATAGTGTCCTAGTAAGGCTTGACAATAACTCTGGTCTTCATCATATTTTGTGATAAAACTTTTGATTTAAAAAAACCTCTGATCTATTTATCATGGCAAATGGATAGAGCTTTCCTGCCTGTTTTCTTTCTTTTCTTTTTTCTTTCTTTCCTTTTTTTTCCTTTGAGCTTAGATTTTTAGAAGCACATATTTAAAAATCAGGTATAAGACTGGATGCAGTGGCTCACGCCTGTAATC)

          chrX:    167,550 out of  3,226,630 (  5.2%) repeat intervals
          chrY:     29,107 out of  3,226,630 (  0.9%) repeat intervals
          chrM:         14 out of  3,226,630 (  0.0%) repeat intervals
   alt contigs:          0 out of  3,226,630 (  0.0%) repeat intervals

Motif size distribution:
          1bp:          0 out of  3,226,630 (  0.0%) repeat intervals
          2bp:    937,030 out of  3,226,630 ( 29.0%) repeat intervals
          3bp:  1,421,462 out of  3,226,630 ( 44.1%) repeat intervals
          4bp:    580,177 out of  3,226,630 ( 18.0%) repeat intervals
          5bp:    175,149 out of  3,226,630 (  5.4%) repeat intervals
          6bp:     55,863 out of  3,226,630 (  1.7%) repeat intervals
       7-24bp:     42,233 out of  3,226,630 (  1.3%) repeat intervals
        25+bp:     14,716 out of  3,226,630 (  0.5%) repeat intervals

Num repeats in reference:
           1x:      7,184 out of  3,226,630 (  0.2%) repeat intervals
           2x:     32,914 out of  3,226,630 (  1.0%) repeat intervals
           3x:  1,776,579 out of  3,226,630 ( 55.1%) repeat intervals
           4x:    620,941 out of  3,226,630 ( 19.2%) repeat intervals
           5x:    331,270 out of  3,226,630 ( 10.3%) repeat intervals
           6x:    138,738 out of  3,226,630 (  4.3%) repeat intervals
           7x:     74,359 out of  3,226,630 (  2.3%) repeat intervals
           8x:     42,620 out of  3,226,630 (  1.3%) repeat intervals
           9x:     31,077 out of  3,226,630 (  1.0%) repeat intervals
       10-15x:    107,308 out of  3,226,630 (  3.3%) repeat intervals
       16-25x:     58,157 out of  3,226,630 (  1.8%) repeat intervals
       26-35x:      5,024 out of  3,226,630 (  0.2%) repeat intervals
       36-50x:        354 out of  3,226,630 (  0.0%) repeat intervals
         51+x:        105 out of  3,226,630 (  0.0%) repeat intervals

Fraction pure bases distribution:
          0.0:        183 out of  3,226,630 (  0.0%) repeat intervals
          0.1:      1,023 out of  3,226,630 (  0.0%) repeat intervals
          0.2:      1,578 out of  3,226,630 (  0.0%) repeat intervals
          0.3:      1,522 out of  3,226,630 (  0.0%) repeat intervals
          0.4:        981 out of  3,226,630 (  0.0%) repeat intervals
          0.5:      4,780 out of  3,226,630 (  0.1%) repeat intervals
          0.6:     10,126 out of  3,226,630 (  0.3%) repeat intervals
          0.7:      3,579 out of  3,226,630 (  0.1%) repeat intervals
          0.8:     20,437 out of  3,226,630 (  0.6%) repeat intervals
          0.9:     36,903 out of  3,226,630 (  1.1%) repeat intervals
          1.0:  3,145,518 out of  3,226,630 ( 97.5%) repeat intervals

Mappability distribution:
          0.0:    102,883 out of  3,226,630 (  3.2%) loci
          0.1:     99,541 out of  3,226,630 (  3.1%) loci
          0.2:    121,854 out of  3,226,630 (  3.8%) loci
          0.3:    175,543 out of  3,226,630 (  5.4%) loci
          0.4:    230,774 out of  3,226,630 (  7.2%) loci
          0.5:    309,700 out of  3,226,630 (  9.6%) loci
          0.6:    227,070 out of  3,226,630 (  7.0%) loci
          0.7:    222,527 out of  3,226,630 (  6.9%) loci
          0.8:    258,323 out of  3,226,630 (  8.0%) loci
          0.9:    512,745 out of  3,226,630 ( 15.9%) loci
          1.0:    965,670 out of  3,226,630 ( 29.9%) loci
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.merge_and_annotate_loci.log)], including:

```
63,170 out of 3,289,806 (  1.9%) adjacent loci that were within 1bp of each other and had the same motif after cyclic shift were merged together

Sources of loci in the final catalog:
   3,220,632 out of 3,289,806 (97.9%) All pure repeats in hg38 (hg38_repeats.motifs_2_to_1000bp.repeats_3x_and_spans_9bp.filtered.json)
      58,447 out of 3,289,806 ( 1.8%) Truth set of TR variants from 51 HPRC assemblies (combined.51_samples.positive_loci.filtered.json)
      10,645 out of 3,289,806 ( 0.3%) Illumina 174k polymorphic loci from 1kGP samples (illumina_variant_catalog.sorted.filtered.json)
          82 out of 3,289,806 ( 0.0%) Known disease-associated loci (variant_catalog_without_offtargets.GRCh38.split.filtered.json)
```


