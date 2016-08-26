#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans'
mpl.rcParams['font.size'] = 11
mpl.rcParams['font.serif'] = 'Calibri'
mpl.rcParams['image.cmap'] = 'coolwarm'
mpl.rcParams['figure.figsize'] = '5,4'

import pylab as pl
from scipy.io import netcdf_file

def plotVar(f,vn):
	I = f.variables[vn]
	x = pl.linspace(1,I.shape[1],I.shape[1])
	y = pl.linspace(1,I.shape[2],I.shape[2])
	
	R = (I.data.min(),I.data.max())

	fig = pl.figure(vn,figsize=[11,5])

	ax = fig.add_subplot(1,2,1,aspect=1.0)
	pc = ax.pcolormesh(x,y,I[0,:,:],vmin=R[0],vmax=R[1])
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel(r'Position $x$ [px]')
	ax.set_ylabel(r'Position $y$ [px]')
	ax.set_title('Image A')

	ax = fig.add_subplot(1,2,2,aspect=1.0)
	pc = ax.pcolormesh(x,y,I[1,:,:],vmin=R[0],vmax=R[1])
	fig.colorbar(pc,ax=ax)
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel(r'Position $x$ [px]')
	ax.set_yticks([])
	ax.set_title('Image B')

	fig.tight_layout(w_pad=-3)

#f = netcdf_file('./results/pixels-3.0/pair-1.nc')
f = netcdf_file('./results/pixels-3.0/vector-1-[1,1|2).nc')
plotVar(f,'I')
plotVar(f,'dudI')
plotVar(f,'dIdU')
plotVar(f,'dvdI')
plotVar(f,'dIdV')
#plotVar(f,'dIdU')
#plotVar(f,'dIdV')
#plotVar(f,'dIdR')
pl.show()
f.close()