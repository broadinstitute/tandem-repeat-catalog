This repo provides steps and tools for defining genome-wide tandem repeat catalogs for use with tools like ExpansionHunter, GangSTR, HipSTR, TRGT, and LongTR. 

---
### Defining a genome-wide tandem repeat catalog

The steps below use examples based on the hg38 human reference, but are applicable to other reference genomes. 

#### Step 1: Detect all perfect (ie. non-interrupted) tandem repeats in the reference genome

We use [colab-repeat-finder](https://github.com/broadinstitute/colab-repeat-finder) - a simple algorithm designed to detect all perfect repeats that satisfy user-specified criteria. 

```
git clone git@github.com:broadinstitute/colab-repeat-finder.git
cd colab-repeat-finder/python
python3 perfect_repeat_finder.py --min-repeats 3 --min-span 9 --min-motif-size 1 --max-motif-size 50  --show-progress-bar /path/to/hg38.fa   --output-prefix perfect_repeats.hg38
```

*NOTE:* We could have, instead, run TandemRepeatFinder (TRF) with very large mismatch and indel penalties to make it return only perfect repeats, but we found that TRF fails to detect ~3% of perfect repeats (2-50bp motifs) for unclear reasons. We therefore use colab-repeat-finder to maximise sensitivity. 


### Step 2: Combine perfect_repeat_finder.py output with any available empirically-defined polymorphic tandem repeat catalogs

While the perfect_repeat_finder.py catalog from step 1 includes all tandem repeat loci that have at least 3 repeats of some motif and span at least 9bp in the reference genome, it misses loci that have 2 or fewer repeats in the reference (while having 3 or more repeats in other genomes within the population). To capture these loci as well, we can merge the catalog from step 1 with catalogs of polymorphic tandem repeat loci that were generated via 
orthogonal methods. 


While doing this, we also fix some representation issues by doing the following:
- collapsing adjacent loci that have the same motif but were reported as two separate loci by colab-repeat-finder due to a single base pair interruption between them.
- filtering out repeats that include non-A,C,G,T characters (such as homopolymer repeats of 'N')


### Step 3: Annotate loci with gene names, etc.

Use the `python3 -u -m str_analysis.annotate_and_filter_str_catalog` script.

### Step 4: Convert the combined catalog into tool-specific catalog formats, while adding tool-specific optimizations to the locus definitions 
