#!/usr/bin/python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

#I/O
full_file = sys.argv[1]
pareto_file = sys.argv[2]

y = np.loadtxt('run/'+full_file)
f = np.loadtxt('run/'+pareto_file)


fig, ax1 = plt.subplots()
#fig.set_size_inches(9, 4)
fig.subplots_adjust(left=0.16,bottom=0.14,right=0.9,top=0.9,wspace=0.5,hspace=0.35)
plt.tick_params(axis='both', labelsize=10, pad=7)
ax1.set_xlabel(r"RMS Bandwidth [%]", size=14)
ax1.set_ylabel(r"Collimated Flux [ph/s]", size=14)
ax1.set_title(r"GA NRB Opt: 200 gen., 50 ind. per gen.", size=14)
ax1.grid(axis='both')

ax1.scatter(y[:,1],-y[:,2],color='b',s=10,label='All GA simulations')
ax1.scatter(f[:,1],-f[:,2],color='r',s=10,label='Pareto Front')

ax1.set_xlim([3,5])
ax1.set_ylim([0,60000])

ax1.legend(loc='upper left',shadow=False,markerscale=1,scatterpoints=1,fontsize=12)
plt.savefig('pareto_Fcol_BW.png',format='png')
plt.show()


fig, ax2 = plt.subplots()
#fig.set_size_inches(9, 4)
fig.subplots_adjust(left=0.16,bottom=0.14,right=0.9,top=0.9,wspace=0.5,hspace=0.35)
plt.tick_params(axis='both', labelsize=10, pad=7)
ax2.set_xlabel(r"$\beta^*_x$ [m]", size=14)
ax2.set_ylabel(r"$\beta^*_y$ [m]", size=14)
ax2.set_title(r"GA NRB Opt: 200 gen., 50 ind. per gen.", size=14)
ax2.grid(axis='both')

ax2.scatter(y[:,4],y[:,5],color='b',s=10,label='All GA simulations')
ax2.scatter(f[:,4],f[:,5],color='r',s=10,label='Pareto Front')

ax2.set_xlim([0.0,20])
ax2.set_ylim([0.0,20])

ax2.legend(loc='upper right',shadow=False,markerscale=1,scatterpoints=1,fontsize=12)
plt.savefig('pareto_betaxy_space.png',format='png')
plt.show()


