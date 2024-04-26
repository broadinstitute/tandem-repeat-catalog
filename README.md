This repo provides steps and tools for defining genome-wide tandem repeat catalogs for use with tools like ExpansionHunter, GangSTR, HipSTR, TRGT, and LongTR. 

---
### Defining a genome-wide tandem repeat catalog

The steps below use examples based on the hg38 human reference, but are applicable to other reference genomes. 

#### Step 1: Detect all perfect (ie. non-interrupted) tandem repeats in the reference genome

We use [colab-repeat-finder](https://github.com/broadinstitute/colab-repeat-finder) - a simple algorithm designed to detect all perfect repeats that satisfy user-specified criteria. 
The reason we use it rather than TandemRepeatFinder (TRF) is that it detects ~3% more repeats (2-50bp motifs) than detected by TRF when TRF is run with very large mismatch and indel penalties that make it return only perfect repeats. 

```
git clone git@github.com:broadinstitute/colab-repeat-finder.git
cd colab-repeat-finder/python
python3 perfect_repeat_finder.py --min-repeats 3 --min-span 9 --min-motif-size 1 --max-motif-size 50  --show-progress-bar /path/to/hg38.fa   --output-prefix perfect_repeats.hg38
```
