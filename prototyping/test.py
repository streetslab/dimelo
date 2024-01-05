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
# plot_enrichment.plot_enrichment_base(['test.bed', 'test.bed'], ['a', 'b'])
# plt.show()
# plt.close()

""" plot_enrichment_profile testing """
# WINDOW_SIZE = 500

# plot_enrichment_profile.plot_enrichment_profile_base(mod_file_names=['test1.fake', 'test2.fake'],
#                                                      bed_file_names=['test1.bed', 'test2.bed'],
#                                                      mod_names=['A', 'C'],
#                                                      window_size=WINDOW_SIZE,
#                                                      sample_names=['sample1', 'sample2'])
# plt.show()
# plt.close()

# plot_enrichment_profile.plot_enrichment_profile_vary_mod(mod_file_name='test.fake',
#                                                          bed_file_name='test.bed',
#                                                          window_size=WINDOW_SIZE,
#                                                          mod_names=['A', 'C'])
# plt.show()
# plt.close()

# plot_enrichment_profile.plot_enrichment_profile_vary_regions(mod_file_name='test.fake',
#                                                              bed_file_names=['test1.bed', 'test2.bed'],
#                                                              mod_name='A',
#                                                              window_size=WINDOW_SIZE,
#                                                              sample_names=['on target', 'off target'])
# plt.show()
# plt.close()

# plot_enrichment_profile.plot_enrichment_profile_vary_experiments(mod_file_names=['test1.fake', 'test2.fake'],
#                                                                  bed_file_name='test.bed',
#                                                                  mod_name='A',
#                                                                  window_size=WINDOW_SIZE,
#                                                                  sample_names=['experiment 1', 'experiment 2'])
# plt.show()
# plt.close()

""" plot_single_reads testing """
plot_single_reads.plot_single_reads_rectangle('test.bed', 'test.bed', ['A', 'C'])
plt.show()
plt.close()

""" random stuff """
# n_reads = 100
# read_len = 500
# x = np.linspace(start=0, stop=1, num=read_len*2)
# y = test_data.expspace_prob(num=read_len*2, a=15)
# z = test_data.expspace_prob(num=read_len*2, a=15, b=0.5)
# plt.plot(x, y)
# plt.plot(x, z)
# plt.show()
# plt.close()