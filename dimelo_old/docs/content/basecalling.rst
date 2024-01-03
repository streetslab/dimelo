Basecalling & Alignment Suggestions
====================================

The dimelo package takes an aligned bam file with modified basecalls as input. Below, we offer some suggestions for performing the upstream basecalling and alignment. 
The best basecalling and alignment methods will vary by use case, and the below suggestions are just what worked well for us.
ONT is also constantly improving basecalling, and these suggestions are likely to become outdated quickly.

Suggestions
------------

1. Model we have used for m6A calling: 

    * `res_dna_r941_min_modbases-all-context_v001 <https://github.com/nanoporetech/rerio/blob/master/basecall_models/res_dna_r941_min_modbases-all-context_v001>`__

2. We generally use Megalodon for both basecalling & alignment. Megalodon requires `Guppy <https://community.nanoporetech.com/downloads>`__. N.B. you may need to downgrade Megalodon and/or Guppy to find a compatible version combination.
    
    * `Megalodon <https://nanoporetech.github.io/megalodon/>`__

3. For exploring repetitve regions of the genome, we've found that Winnowmap performs better than the aligner Megalodon uses, `minimap2 <https://github.com/lh3/minimap2>`__. For repetitive regions, we use Guppy for basecalling and Winnowmap for alignment.
    
    * `Guppy <https://community.nanoporetech.com/downloads>`__
    * `Winnowmap <https://github.com/marbl/Winnowmap>`__

4. If using Guppy & Winnowmap, the resulting bam files must be combined to create a single bam with the modified basecalls from Guppy and the mapping information from Winnowmap. 
   
    * `hybrid_guppy_winnnowmap_bam_creation <https://github.com/amaslan/dimelo-seq/blob/main/ctcf_and_h3k9me3/hybrid_guppy_winnnowmap_bam_creation>`__

5. We perform basecalling separately from the sequencing run and use an EC2 instance (g4dn.metal) with multiple GPUs to speed up basecalling significantly. 