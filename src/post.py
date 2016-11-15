#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans'
mpl.rcParams['font.size'] = 11
mpl.rcParams['font.serif'] = 'Calibri'
mpl.rcParams['image.cmap'] = 'coolwarm'
mpl.rcParams['figure.figsize'] = '4,4'

import pylab as pl
from scipy.io import netcdf_file

Nsets = 7

seq = pl.get_cmap('Greens')
div = pl.get_cmap('PRGn')

capkeys = []
captions = {}


Ncf = 30

#===================#
#= Plot Image Pair =#
#===================#

def doImage(I,nm,sn):
	x = pl.linspace(1,I.shape[0],I.shape[0])
	y = pl.linspace(1,I.shape[1],I.shape[1])
	
	fig = pl.figure(figsize=(4,4),tight_layout=True)
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	
	ax.pcolormesh(x,y,I[:,:],cmap=seq)
	#ax.contourf(x,y,I[:,:],Ncf,cmap=seq)
	
	ax.axhline(47.5,color='k',ls='--',lw=2)
	ax.axhline(96.5,color='k',ls='--',lw=2)
	ax.axvline(47.5,color='k',ls='--',lw=2)
	ax.axvline(96.5,color='k',ls='--',lw=2)
	
	ax.plot([72.5],[72.5],'kx',ms=10,mew=3)
	
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel('Position $x$ [px]')
	ax.set_ylabel('Position $y$ [px]')
	
	fn = '%s-%s.pdf'%(sn,nm)
	fig.savefig(fn)
	capkeys.append(fn)
	captions[fn] = 'Image Intensity $I$ $[%s]$'%nm
	pl.close(fig)

def doPair(fn,sn):
	f = netcdf_file(fn)
	Ia = f.variables['I'][0,:,:].copy()
	Ib = f.variables['I'][1,:,:].copy()
	f.close()
	
	doImage(Ia,'A',sn)
	doImage(Ib,'B',sn)

#===============#
#= Plot Vector =#
#===============#

def doVectorVar(V,vn,fvn,nm,sn,sym=True):
	x = pl.linspace(48,48+V.shape[0]-1,V.shape[0])
	y = pl.linspace(48,48+V.shape[1]-1,V.shape[1])
	
	fig = pl.figure(figsize=(4,4),tight_layout=True)
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	
	if sym:
		cm = div
		M = pl.median(V)
		R = abs(V-M).max()
		vr = [M-R,M+R]
	else:
		cm = seq
		vr = [V.min(),V.max()]
	
	ax.pcolormesh(x,y,V[:,:],cmap=cm,vmin=vr[0],vmax=vr[1])
	#ax.contourf(x,y,V[:,:],Ncf,cmap=cm,vmin=vr[0],vmax=vr[1])
	
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel('Position $x$ [px]')
	ax.set_ylabel('Position $y$ [px]')
	
	fn = '%s-%s-%s.pdf'%(sn,vn,nm)
	fig.savefig(fn)
	capkeys.append(fn)
	captions[fn] = 'Vector\'s Regional %s $[%s]$'%(fvn,nm)
	pl.close(fig)

def doVector(fn,sn):
	vars = ['I','dIdR','dIdN','dIdU','dIdUx','dIdUy','dIdV','dIdVx','dIdVy','dudI','dvdI']
	fancyVars = {}
	fancyVars['I'] = '$I$'
	fancyVars['dIdU']  = r'$\frac{\partial I}{\partial U}$'
	fancyVars['dIdUx'] = r'$\frac{\partial I}{\partial U_x}$'
	fancyVars['dIdUy'] = r'$\frac{\partial I}{\partial U_y}$'
	fancyVars['dIdV']  = r'$\frac{\partial I}{\partial V}$'
	fancyVars['dIdVx'] = r'$\frac{\partial I}{\partial V_x}$'
	fancyVars['dIdVy'] = r'$\frac{\partial I}{\partial V_y}$'
	fancyVars['dIdR']  = r'$\frac{\partial I}{\partial R}$'
	fancyVars['dIdN']  = r'$\frac{\partial I}{\partial N}$'
	fancyVars['dudI']  = r'$\frac{\partial u}{\partial I}$'
	fancyVars['dvdI']  = r'$\frac{\partial v}{\partial I}$'
	
	f = netcdf_file(fn)
	N = f.variables['I'][0,:,:].shape
	
	A = {}
	B = {}
	for v in vars:
		A[v] = f.variables[v][0,:,:].copy()
		B[v] = f.variables[v][1,:,:].copy()
	
	f.close()
	
	for v in vars[:3]:
		fv = fancyVars[v]
		doVectorVar(A[v],v,fv,'A',sn,sym=False)
		doVectorVar(B[v],v,fv,'B',sn,sym=False)
	
	for v in vars[3:]:
		fv = fancyVars[v]
		doVectorVar(A[v],v,fv,'A',sn)
		doVectorVar(B[v],v,fv,'B',sn)

#========================#
#= Plot Correlation Map =#
#========================#

def doMapVar(V,vn,fvn,sn,sym=True):
	xe = (V.shape[0]-1)/2
	ye = (V.shape[1]-1)/2
	
	x = pl.linspace(-xe,xe,V.shape[0])
	y = pl.linspace(-ye,ye,V.shape[1])
	
	fig = pl.figure(figsize=(4,4.25),tight_layout=True)
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	
	if sym:
		cm = div
		M = pl.median(V)
		R = abs(V-M).max()
		vr = [M-R,M+R]
	else:
		cm = seq
		vr = [V.min(),V.max()]
	
	ax.pcolormesh(x,y,V[:,:],cmap=cm,vmin=vr[0],vmax=vr[1])
	#ax.contourf(x,y,V[:,:],Ncf,cmap=cm,vmin=vr[0],vmax=vr[1])
	
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel('Position $x$ [px]')
	ax.set_ylabel('Position $y$ [px]')
	
	fn = '%s-%s.pdf'%(sn,vn)
	fig.savefig(fn)
	capkeys.append(fn)
	captions[fn] = 'Correlation Map %s'%(fvn)
	pl.close(fig)

def doMap(fn,sn):
	vars = ['I','dIdR','dIdN','dIdU','dIdUx','dIdUy','dIdV','dIdVx','dIdVy']
	fancyVars = {}
	fancyVars['I'] = '$C$'
	fancyVars['dIdU']  = r'$\frac{\partial C}{\partial U}$'
	fancyVars['dIdUx'] = r'$\frac{\partial C}{\partial U_x}$'
	fancyVars['dIdUy'] = r'$\frac{\partial C}{\partial U_y}$'
	fancyVars['dIdV']  = r'$\frac{\partial C}{\partial V}$'
	fancyVars['dIdVx'] = r'$\frac{\partial C}{\partial V_x}$'
	fancyVars['dIdVy'] = r'$\frac{\partial C}{\partial V_y}$'
	fancyVars['dIdR']  = r'$\frac{\partial C}{\partial R}$'
	fancyVars['dIdN']  = r'$\frac{\partial C}{\partial N}$'
	
	f = netcdf_file(fn)
	N = f.variables['I'][:,:].shape
	
	A = {}
	for v in vars:
		A[v] = f.variables[v][0,:,:].copy()
	
	f.close()
	
	for v in vars[:3]:
		fv = fancyVars[v]
		doMapVar(A[v],v,fv,sn,sym=False)
	
	for v in vars[3:]:
		fv = fancyVars[v]
		doMapVar(A[v],v,fv,sn)

