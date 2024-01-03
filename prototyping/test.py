import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import utils
import test_data
import plot_enrichment
import plot_enrichment_profile
import plot_single_reads

""" plot_enrichment testing """
plot_enrichment.plot_enrichment_base(['test.bed', 'test.bed'], ['a', 'b'])
plt.show()
plt.close()

""" plot_enrichment_profile testing """
plot_enrichment_profile.plot_enrichment_profile_base(mod_file_names=['test1.bed', 'test2.bed'],
                                                     bed_file_names=['test1.bed', 'test2.bed'],
                                                     mod_names=['A', 'C'],
                                                     sample_names=['sample1', 'sample2'])
plt.show()
plt.close()
plot_enrichment_profile.plot_enrichment_profile_vary_mod(mod_file_name='test.bed',
                                                         bed_file_name='test.bed',
                                                         mod_names=['A', 'C'])
plt.show()
plt.close()
plot_enrichment_profile.plot_enrichment_profile_vary_regions(mod_file_name='test.bed',
                                                             bed_file_names=['test1.bed', 'test2.bed'],
                                                             mod_name='A',
                                                             sample_names=['on target', 'off target'])
plt.show()
plt.close()

plot_enrichment_profile.plot_enrichment_profile_vary_experiments(mod_file_names=['test1.bed', 'test2.bed'],
                                                                 bed_file_name='test.bed',
                                                                 mod_name='A',
                                                                 sample_names=['experiment 1', 'experiment 2'])
plt.show()
plt.close()

""" plot_single_reads testing """
plot_single_reads.plot_single_reads_rectangle('test.bed', 'test.bed', ['A', 'C'])
plt.show()
plt.close()