import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from . import utils
from . import test_data
from . import plot_enrichment
from . import plot_enrichment_profile
from . import plot_single_reads

""" plot_enrichment testing """
plot_enrichment.plot_enrichment_base(mod_file_names=['test.fake', 'test.fake'],
                                     bed_file_names=['test.bed', 'test.bed'],
                                     mod_names=['A', 'A'],
                                     sample_names=['a', 'b'])
plt.show()
plt.close()

""" plot_enrichment_profile testing """
WINDOW_SIZE = 500
SMOOTH_WINDOW=50

plot_enrichment_profile.plot_enrichment_profile_base(mod_file_names=['test1.fake', 'test2.fake'],
                                                     bed_file_names=['test1.bed', 'test2.bed'],
                                                     mod_names=['A', 'C'],
                                                     window_size=WINDOW_SIZE,
                                                     sample_names=['sample1', 'sample2'],
                                                     smooth_window=SMOOTH_WINDOW)
plt.show()
plt.close()

plot_enrichment_profile.plot_enrichment_profile_vary_mod(mod_file_name='test.fake',
                                                         bed_file_name='test.bed',
                                                         window_size=WINDOW_SIZE,
                                                         mod_names=['A', 'C'],
                                                         smooth_window=SMOOTH_WINDOW)
plt.show()
plt.close()

plot_enrichment_profile.plot_enrichment_profile_vary_regions(mod_file_name='test.fake',
                                                             bed_file_names=['test1.bed', 'test2.bed'],
                                                             mod_name='A',
                                                             window_size=WINDOW_SIZE,
                                                             sample_names=['on target', 'off target'],
                                                             smooth_window=SMOOTH_WINDOW)
plt.show()
plt.close()

plot_enrichment_profile.plot_enrichment_profile_vary_experiments(mod_file_names=['test1.fake', 'test2.fake'],
                                                                 bed_file_name='test.bed',
                                                                 mod_name='A',
                                                                 window_size=WINDOW_SIZE,
                                                                 sample_names=['experiment 1', 'experiment 2'],
                                                                 smooth_window=SMOOTH_WINDOW)
plt.show()
plt.close()

""" plot_single_reads testing """
# plot_single_reads.plot_single_reads_rectangle('test.bed', 'test.bed', ['A', 'C'])
# plt.show()
# plt.close()

""" random stuff """
