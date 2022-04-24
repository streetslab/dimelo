r"""
===============
dimelo module
===============
.. currentmodule:: dimelo

dimelo allows you to perform quality control and plot modified bases from bam files.

.. automodule:: parse_bam
.. automodule:: plot_browser
.. automodule:: plot_enrichment
.. automodule:: plot_enrichment_profile
.. automodule:: qc_report

"""
from .parse_bam import parse_bam
from .plot_browser import plot_browser
from .plot_enrichment import plot_enrichment
from .plot_enrichment_profile import plot_enrichment_profile
from .qc import qc_report
