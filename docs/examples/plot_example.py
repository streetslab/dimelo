"""
"This" is my example-script
===========================

This example doesn't do much, it just makes a simple plot
"""
# %%
# This is a section header
# ------------------------
# This is the first section!
# The `#%%` signifies to Sphinx-Gallery that this text should be rendered as
# rST and if using one of the above IDE/plugin's, also signifies the start of a
# 'code block'.

# This line won't be rendered as rST because there's a space after the last block.
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel(r'$x$')
plt.ylabel(r'$\sin(x)$')
# To avoid matplotlib text output
plt.show()
# This is the end of the 'code block' (if using an above IDE). All code within
# this block can be easily executed all at once.

# %%
# This is another section header
# ------------------------------
#
# In the built documentation, it will be rendered as rST after the code above!
# This is also another code block.
print("here")