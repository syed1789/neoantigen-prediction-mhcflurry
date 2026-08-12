## Difference between v1.2 and v1.2.1 VCF files

The two files `high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz` and `high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz` were filtered for calls in the consensus callable regions 
(i.e., `mergeBed.MajorityAlignersCallable.bed` as described [here](https://sites.google.com/view/seqc2/home/identify-high-confidence-regions).

However, the high-confidence regions (i.e., `High-Confidence_Regions_v1.2.bed`) subtracted additional 20-bp regions spanning the low-confidence and likely false positive positions in the all-inclusive call set.
There are 113 SNVs and 297 indels in those subtracted regions, and they were left in `high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz` and `high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz`.

We used bedtools intersect to create `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` and `high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz`. All the calls there are inside `High-Confidence_Regions_v1.2.bed`.
