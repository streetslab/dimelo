# Test files
## Both bams are free pA-Hia5 control in GM12878 cells
- Example bam from megalodon: mod_mappings_subset.bam
- Example hybrid bam from winnowmap & guppy merge: winnowmap_guppy_merge_subset.bam


## Test files created using samtools subsample
File size limit for github is 100 MB and <50 MB is recommended:
```
samtools view -s 0.01 -@ 10 prod_free_Hia5_mod_mappings.sorted.bam -o mod_mappings_subset.bam

samtools view -s 0.004 -@ 10 20210512_6_winnnowmap_guppy_merge.sorted.bam -o winnowmap_guppy_merge_subset.bam
```