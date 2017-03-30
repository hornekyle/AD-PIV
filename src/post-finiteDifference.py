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

cm = pl.get_cmap('PRGn')

def getIntensity(fn):
	fh = netcdf_file(fn)
	O = {}
	N = fh.variables['I'].shape[2:]
	for v in ['I','dIdU']:
		O[v] = pl.empty(N)
		O[v][:,:] = fh.variables[v][0,0,:,:]
	fh.close()
	return O

def getVector(fn):
	fh = netcdf_file(fn)
	O = {}
	N = fh.variables['u'].shape[1:]
	for v in ['u','dudU','v','dvdU']:
		O[v] = pl.empty(N)
		O[v][:,:] = fh.variables[v][0,:,:]
	fh.close()
	return O

def getSeries(disps,nfn,getter):
	O = {}
	for d in disps:
		off = int(round(10.0*(d-5.0)))
		fn = './results/finiteDifference/%2.1f/%s'%(d,nfn)
		O[off] = getter(fn)
	return O

def plotField(F,T):
	fig = pl.figure(tight_layout=True,figsize=(5,4))
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	im = ax.imshow(F,cmap=cm,interpolation='nearest')
	ax.set_xlabel('x-coordinate [px]')
	ax.set_ylabel('y-coordinate [px]')
	fig.colorbar(im,ax=ax,label=T)
	return fig

def rootSumSquares(F):
	return pl.sqrt((F**2).sum())

def plotConvergence(ax,h0,e0,c,s,t):
	h = h0[:]/h0[-1]
	e = e0[:]/e0[-1]
	m,b = pl.polyfit(pl.log(h),pl.log(e),1)
	ax.loglog(h,pl.exp(b)*h**m,'%s%s'%(c,s[0]),label=r'$\epsilon \propto O(h^{%.2f})$'%m,zorder=1)
	ax.loglog(h,e,'k%s'%s[1],label='%s'%t,zorder=2)

def computeConvergence(S,v,d):
	D_AD = S[0][d]
	D_FD = {}
	offsets = range(1,int((disps.size-1)/2+1))
	for k in offsets: D_FD[k]  = (-S[-k][v]+S[k][v])/(2*k*h0)
	spacing = pl.array([ h0*k for k in offsets])
	errors  = pl.array([ rootSumSquares(D_FD[k]-D_AD) for k in offsets ])
	return spacing,errors

disps = pl.linspace(4.4,5.6,13)
h0 = disps[1]-disps[0]
Su = getSeries(disps,'vectors.nc',getVector)
SI = getSeries(disps,'vector-[1,1|1].nc',getIntensity)

hu,eu = computeConvergence(Su,'u','dudU')
hI,eI = computeConvergence(SI,'I','dIdU')

fig = pl.figure(figsize=(4,3),tight_layout=True)
ax = fig.add_subplot(1,1,1)
plotConvergence(ax,hI,eI,'C0',['--','s'],r'$\partial I / \partial U$')
plotConvergence(ax,hu[:-2],eu[:-2],'C1',[':','o'],r'$\partial u / \partial U$')
ax.legend(loc='upper center',ncol=2)
ax.set_xlim(10**(-1.0),10**(0.2))
ax.set_ylim(10**(-1.8),10**(0.8))
ax.set_xlabel('Normalized Spacing $h$ [-]')
ax.set_ylabel('Normalized Error $\\epsilon$ [-]')

fig.savefig('./figures/finiteDifference/convergence.svg')
pl.close(fig)
