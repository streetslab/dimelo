dimelo
========

a python package & command line tool for analyzing DiMeLo-seq data

Functions
--------------
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| function                                        | type                  | input files                         | main output file(s)                                                                                                            |
+=================================================+=======================+=====================================+============================================+===================================================================================+
|:doc:`content/reference/qc_report`               | QC                    | 1 or more BAM                       | summary report with figures & table for read length, mapping quality, basecall quality, etc.                                   |
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`content/reference/plot_browser`            | visualization         | 1 or more BAM                       | single-molecule browser as html or pdf; aggregate browser plots and coverage                                                   |
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`content/reference/plot_enrichment`         | visualization         | 1 or more BAM, 1 or more BED        | barplot of fraction modified bases within region(s) defined by BED(s)                                                          |
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`content/reference/plot_enrichment_profile` | visualization         | 1 or more BAM, 1 or more BED        | profile, single molecule plots, base abundance; option to overlay if input multiple BAM/BED                                    |
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`content/reference/parse_bam`               | data formatting       | 1 BAM                               | sql database (to allow for custom visualization)                                                                               |
+-------------------------------------------------+-----------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+

Resources
-----------------------

:doc:`content/installation` 

:doc:`auto_examples/index`

:doc:`content/basecalling` 

:doc:`content/tutorial` 

`Open Source Code <https://github.com/amaslan/dimelo>`__

`DiMeLo-seq Protocol <https://www.protocols.io/view/dimelo-seq-directed-methylation-with-long-read-seq-n2bvjxe4wlk5/v2/materials>`__


.. toctree::
   :maxdepth: 2
   :hidden:
   
   content/installation
   content/reference/tools
   auto_examples/index
   content/basecalling
   content/tutorial


