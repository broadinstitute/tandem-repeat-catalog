

### Repeat Catalog

This repo provides a general purpose genome-wide TR catalog for genotyping TRs in short read or long read sequencing data. 

[Release v0.9](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v0.9) is now available on the [[Releases]](https://github.com/broadinstitute/tandem-repeat-catalogs/releases) page and includes the following files:

<br />
<table>
<tr><td><b>File</b></td><td><b>Description</b></td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.merged_and_annotated.json.gz</td><td>Includes all annotations, and can be used directly with <a href="https://github.com/Illumina/ExpansionHunter">ExpansionHunter</a></td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.bed.gz</td><td>Sorted and indexed BED file for IGV</td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.TRGT.bed</td><td>For use with <a href="https://github.com/PacificBiosciences/trgt">TRGT</a></td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.LongTR.bed</td><td>For use with <a href="https://github.com/gymrek-lab/LongTR">LongTR</a></td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.GangSTR.bed</td><td>For use with <a href="https://github.com/gymreklab/GangSTR">GangSTR</a></td></tr>
<tr><td>repeat_catalog.3x_and_9bp.2_to_1000bp_motifs.hg38.HipSTR.bed</td><td>For use with <a href="https://github.com/HipSTR-Tool/HipSTR">HipSTR</a></td></tr>
</table>
<br />

### Goals

- Create a catalog that is as sensitive (includes all polymorphic TR loci in the human genome) and specific (excludes non-polymorphic loci) as possible
- Provide rich annotations
- Share the catalog in formats that can be used directly with widely-used TR genotyping tools

Although our initial focus is on the human genome, we'd also love to extend this work to plants and other species. Please consider creating a GitHub issue or reaching out by email if you are interested in this.

### Background

Tandem repeats (TRs) are regions of the genome that consist of consecutive copies of some motif sequence. For example, `CAGCAGCAG` is a tandem repeat of the `CAG` motif. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs, which specify the genomic start and end coordinates of each tandem repeat region, as well as the one or more motifs that repeat there. 

For example, if a hypothetical region at the beginning of `chrX` had the following nucleotide sequence:  
`ATCAGTAGA ATATATATAT CAGACAGCAGCAG TGAGTGCGTAC...`  
it could be represented in a repeat catalog as two entries:  
`chrX:10-19 (AT)*`  
`chrX:20-32 (CAG)*`   
indicating that a repeat of the `AT` motif occurs between positions 10 and 19 (inclusive), and of the `CAG` motif between positions 20 and 32.

