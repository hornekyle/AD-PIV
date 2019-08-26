#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11
mpl.rcParams['font.serif'] = 'palatino'
mpl.rcParams['font.sans-serif'] = 'avant guard'
mpl.rcParams['text.usetex'] = 'yes'
mpl.rcParams['image.cmap'] = 'viridis'

import pylab as pl
from scipy.io import netcdf_file

def plotPair(fn):
	fh = netcdf_file(fn)
	
	Nx = fh.variables['x'].shape[0]
	x  = pl.linspace(1,Nx,Nx)
	
	Ny = fh.variables['y'].shape[0]
	y  = pl.linspace(1,Ny,Ny)
	
	I = pl.empty( (Nx,Ny) )
	I[:,:] = fh.variables['I'][0,1,:,:]
	
	fig = pl.figure(figsize=(4,3),tight_layout=True)
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	cf = ax.pcolorfast(x,y,I,cmap='BuPu')
	
	for k in [48,96]:
		ax.axhline(k,color='k',linestyle='--')
		ax.axvline(k,color='k',linestyle='--')
	
	ax.set_xticks([1,24,48,72,96,120,144])
	ax.set_yticks([1,24,48,72,96,120,144])
	ax.set_xlabel('Position $x$ $/$ px')
	ax.set_ylabel('Position $y$ $/$ px')
	fig.colorbar(cf,ax=ax,label='Intensity $I$ $/$ 1')
	fig.savefig('./figures/example/pair.svg')
	
	fh.close()

def plotVector(fn):
	fh = netcdf_file(fn)
	
	Nx = fh.variables['x'].shape[0]
	x  = pl.linspace(1,Nx,Nx)+47.5
	
	Ny = fh.variables['y'].shape[0]
	y  = pl.linspace(1,Ny,Ny)+47.5
	
	v = pl.empty( (Nx,Ny) )
	
	vns = ['I','dIdU','dIdV','dIdR','dudI','dvdI']
	vds = ['Intensity $I$ [-]',
	       'Intensity Derivative $\partial_U (I)$ $/$ px$^{-1}$',
	       'Intensity Derivative $\partial_V (I)$ $/$ px$^{-1}$',
	       'Intensity Derivative $\partial_R (I)$ $/$ px$^{-1}$',
	       'Velocity Derivative $\partial_I (u)$ $/$ px',
	       'Velocity Derivative $\partial_I (v)$ $/$ px'
	       ]
	vDs = [False,True,True,False,True,True]
	for k in range(len(vns)):
		v[:,:] = fh.variables[ vns[k] ][0,1,:,:]
		
		fig = pl.figure(figsize=(4,3),tight_layout=True)
		ax = fig.add_subplot(1,1,1,aspect=1.0)
		if vDs[k]:
			ex = max(abs(v.min()),abs(v.max()))
			cf = ax.pcolorfast(x,y,v,cmap='PuOr',vmin=-ex,vmax=ex)
		else:
			cf = ax.pcolorfast(x,y,v,cmap='BuPu')
		
		ax.set_xticks([48,72,96])
		ax.set_yticks([48,72,96])
		ax.set_xlabel('Position $x$ $/$ px')
		ax.set_ylabel('Position $y$ $/$ px')
		fig.colorbar(cf,ax=ax,label='%s'%vds[k])
		fig.savefig('./figures/example/vector-%s.svg'%vns[k])
	
	fh.close()

def plotMap(fn):
	fh = netcdf_file(fn)
	
	Nx = fh.variables['x'].shape[0]
	x  = pl.linspace(-(Nx-1)/2,(Nx-1)/2,Nx)
	
	Ny = fh.variables['y'].shape[0]
	y  = pl.linspace(-(Ny-1)/2,(Ny-1)/2,Ny)
	
	v = pl.empty( (Nx,Ny) )
	
	vns = ['I','dIdU','dIdV','dIdR']
	vds = ['Correlation $C$ $/$ 1',
	       'Correlation Derivative $\partial_U (C)$ $/$ px$^{-1}$',
	       'Correlation Derivative $\partial_V (C)$ $/$ px$^{-1}$',
	       'Correlation Derivative $\partial_R (C)$ $/$ px$^{-1}$'
	       ]
	vDs = [False,True,True,False]
	for k in range(len(vns)):
		v[:,:] = fh.variables[ vns[k] ][0,:,:]
		
		fig = pl.figure(figsize=(4,3),tight_layout=True)
		ax = fig.add_subplot(1,1,1,aspect=1.0)
		if vDs[k]:
			ex = max(abs(v.min()),abs(v.max()))
			cf = ax.pcolorfast(x,y,v,cmap='PuOr',vmin=-ex,vmax=ex)
		else:
			cf = ax.pcolorfast(x,y,v,cmap='BuPu')
		
		ax.set_xlabel('Position $x$ $/$ px')
		ax.set_ylabel('Position $y$ $/$ px')
		fig.colorbar(cf,ax=ax,label='%s'%vds[k])
		fig.savefig('./figures/example/map-%s.svg'%vns[k])
	
	fh.close()

plotPair('./results/example/pair.nc')
plotVector('./results/example/vector-[1,1|2].nc')
plotMap('./results/example/map-[1,1|1].nc')
