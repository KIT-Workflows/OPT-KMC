# -*- coding: utf-8 -*-
"""analysis_opt-PSO.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1TzZ_5ICH4YCvL_s4UjXOrustPcRHcYNw

### Libraries
"""

import os,sys
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

plt.style.use("data_input/figure_style.mplstyle")

"""### Upload data"""

loss_func = pd.read_csv('analysis_opt-PSO/Loss_function.csv', header=None)
en_bar = pd.read_csv('analysis_opt-PSO/Energy_barrier.csv', header=None)

"""### Plot results"""

fig1 = plt.figure(1, figsize = (4.5,4.5), dpi=100)
ax1 = fig1.add_subplot(111)

ax1.plot(loss_func.iloc[:, 0], color='#088C9C')

ax1.set_yscale('log')
ax1.set_ylabel('Loss function')
ax1.set_xlabel('# itarations')

# fig1.savefig('analysis_opt-PSO/loss_function.png', dpi=200)

# name_par = ['EC-Li', 'P', 'I', 'F', 'O', 'diff_O','diff_P', 'diff_EC-Li', 'diff_Li', 'conc_Li']

fig2 = plt.figure(1, figsize=(7.5,3.5), dpi=100)
grid = fig2.add_gridspec(1, 2, hspace=0.2, wspace=0.3)
ax1 = fig2.add_subplot(grid[0])
ax2 = fig2.add_subplot(grid[1])


im = ax1.scatter(en_bar.iloc[:, 0], en_bar.iloc[:, 7], s=5, c=np.log(loss_func[0]))
im = ax2.scatter(en_bar.iloc[:, 1], en_bar.iloc[:, 6], s=5, c=np.log(loss_func[0]))


ax1.set_xlabel('$E_{a}$ [eV] Reaction- $cEC^{-}$')
ax1.set_ylabel('$E_{a}$ [eV] Diffusion- $cEC^{-}$')

ax2.set_xlabel('$E_{a}$ [eV] Reaction- $oEC^{-}$')
ax2.set_ylabel('$E_{a}$ [eV] Diffusion- $oEC^{-}$')

clb = fig2.colorbar(im, orientation='vertical')
clb.set_label('Loss function')

fig3 = plt.figure(1, figsize=(7.5,3.5), dpi=100)
grid = fig3.add_gridspec(1, 2, hspace=0.2, wspace=0.3)
ax1 = fig3.add_subplot(grid[0])
ax2 = fig3.add_subplot(grid[1])



im = ax1.scatter(en_bar.iloc[:, 2], en_bar.iloc[:, 3], s=5, c=np.log(loss_func[0]))
im = ax2.scatter(en_bar.iloc[:, 4], en_bar.iloc[:, 5], s=5, c=np.log(loss_func[0]))


ax1.set_xlabel('$E_{a}$ [eV] Reaction- $CO3^{2-}-Intermediate$')
ax1.set_ylabel('$E_{a}$ [eV] Reaction- $CO3^{2-}$')

ax2.set_xlabel('$E_{a}$ [eV] Reaction- $EDC^{2-}$')
ax2.set_ylabel('$E_{a}$ [eV] Diffusion- $EDC^{2-}$')

clb = fig3.colorbar(im, orientation='vertical')
clb.set_label('Loss function')

fig4 = plt.figure(1, figsize=(3.5,3.5), dpi=100)
grid = fig4.add_gridspec(1, 1, hspace=0.2, wspace=0.3)
ax1 = fig4.add_subplot(grid[0])


im = ax1.scatter(en_bar.iloc[:, 9], en_bar.iloc[:, 8], s=5, c=np.log(loss_func[0]))

ax1.set_xlabel('$Li^{+}$ Concentration')
ax1.set_ylabel('$E_{a}$ [eV] Diffusion- $Li^{+}$')


clb = fig4.colorbar(im, orientation='vertical')
clb.set_label('Loss function')

# fig2.savefig('analysis_opt-PSO/energy_barrier.png',dpi=200)

