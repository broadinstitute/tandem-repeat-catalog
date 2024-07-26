

### Repeat Catalog

This repo provides a general purpose genome-wide TR catalog for use in TR genotyping studies involving short read or long read sequencing data. 

**Releases**

[v0.9](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v0.9) is now available on the [[Releases]](https://github.com/broadinstitute/tandem-repeat-catalogs/releases) page

### Goals

- Provide a catalog that is as sensitive (includes all polymorphic TR loci in the human genome) and specific (excludes non-polymorphic loci) as possible
- Provide rich annotations
- Determine best practices for how to define TR loci for specific TR genotyping tools, and then share the catalog in formats that can be used directly with:
  * [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)
  * [TRGT](https://github.com/PacificBiosciences/trgt)
  * [LongTR](https://github.com/gymrek-lab/LongTR)
  * [HipSTR](https://github.com/HipSTR-Tool/HipSTR)
  * [GangSTR](https://github.com/gymreklab/GangSTR)
  * and other TR genotyping tools

Although our initial focus is on the human genome, we'd also love to extend this work to plants and other species. Please consider creating a GitHub issue or reaching out by email if you are interested in this.

### Authors

* Ben Weisburd, Broad Institute
* Egor Dolzhenko, PacBio

### Background

Tandem repeats (TRs) are regions of the genome that consist of consecutive copies of some motif sequence. For example, `CAGCAGCAG` is a tandem repeat of the `CAG` motif. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs, which specify the genomic start and end coordinates of each tandem repeat region, as well as the one or more motifs that repeat there. 

For example, if a hypothetical region at the beginning of `chrX` had the following nucleotide sequence:  
`ATCAGTAGA ATATATATAT CAGACAGCAGCAG TGAGTGCGTAC...`  
it could be represented in a repeat catalog as two entries:  
`chrX:10-19 (AT)*`  
`chrX:20-32 (CAG)*`   
indicating that a repeat of the `AT` motif occurs between positions 10 and 19 (inclusive), and of the `CAG` motif between positions 20 and 32.

