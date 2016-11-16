#!/usr/bin/env python

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans'
mpl.rcParams['font.size'] = 11
mpl.rcParams['font.serif'] = 'Calibri'
mpl.rcParams['image.cmap'] = 'viridis'
mpl.rcParams['figure.figsize'] = '5,4'

import pylab as pl
from scipy.io import netcdf_file

seq = pl.get_cmap('viridis')
div = pl.get_cmap('coolwarm')

Ncf = 20+1
Np  = 128

def plotImage(I,nm,sn,vr=[0,1]):
	x = pl.linspace(1,I.shape[0],I.shape[0])
	y = pl.linspace(1,I.shape[1],I.shape[1])
	
	fig = pl.figure(figsize=(5,4),tight_layout=True)
	ax = fig.add_subplot(1,1,1,aspect=1.0)
	
	c = ax.pcolormesh(x,y,I[:,:],cmap=seq,vmin=vr[0],vmax=vr[1])
	#c = ax.contourf(x,y,I[:,:],cmap=seq,levels=pl.linspace(vr[0],vr[1],Ncf))
	
	ax.axhline(47.5,color='k',ls='--',lw=2)
	ax.axhline(96.5,color='k',ls='--',lw=2)
	ax.axvline(47.5,color='k',ls='--',lw=2)
	ax.axvline(96.5,color='k',ls='--',lw=2)
	
	ax.plot([72.5],[72.5],'kx',ms=10,mew=3)
	
	ax.set_xlim(x.min(),x.max())
	ax.set_ylim(y.min(),y.max())
	ax.set_xlabel('Position $x$ [px]')
	ax.set_ylabel('Position $y$ [px]')
	
	fig.colorbar(c,ax=ax)
	
	fn = '%s-%s.pdf'%(sn,nm)
	fig.savefig(fn)
	pl.close(fig)

def doPair(fns):
	f = netcdf_file(fns[0])
	Ia = f.variables['I'][0,:,:].copy()
	Ib = f.variables['I'][1,:,:].copy()
	f.close()
	
	IA = pl.empty( [Ia.shape[0] , Ia.shape[1] , len(fns) ] )
	IB = pl.empty( [Ia.shape[0] , Ia.shape[1] , len(fns) ] )
	
	for k in range( len(fns) ):
		fn = fns[k]
		f = netcdf_file(fn)
		IA[:,:,k] = f.variables['I'][0,:,:].copy()
		IB[:,:,k] = f.variables['I'][1,:,:].copy()
		f.close()
	
	plotImage(IA.mean(2),'A','mean')
	plotImage(IB.mean(2),'B','mean')
	
	plotImage(IA.std(2),'A','std')
	plotImage(IB.std(2),'B','std')

def plotHistogram(I,name,title='',Nb=None):
	if Nb==None:
		Nb = int(pl.sqrt(I.size))
	chance,edges = pl.histogram(I,Nb,density=True)
	edgesLeft  = edges[:-1]
	edgesRight = edges[1:]
	centers = (edgesLeft+edgesRight)/2.0
	widths = edgesRight-edgesLeft
	
	smap = mpl.cm.ScalarMappable()
	smap.set_array(chance)
	colors = smap.to_rgba(chance)
	
	fig = pl.figure(tight_layout=True,figsize=(5,4))
	ax = fig.add_subplot(1,1,1)
	ax.bar(centers,chance,0.8*widths,lw=2,align='center',color=colors)
	ax.set_xlabel(title)
	ax.set_yticks([])
	fig.savefig('%s.pdf'%name)
	pl.close(fig)

def doVector(fns):
	Nf = len(fns)
	Np = 4
	
	keys = ['u','dudU','dudV','dudUx','dudUy','dudVx','dudVy','v','dvdU','dvdV','dvdUx','dvdUy','dvdVx','dvdVy']
	titles = {}
	titles['u']     = 'Displacement $u$ [px]'
	titles['dudU']  = 'Displacement Derivative $\\partial u / \\partial U$ [px/px]'
	titles['dudV']  = 'Displacement Derivative $\\partial u / \\partial V$ [px/px]'
	titles['dudUx'] = 'Displacement Derivative $\\partial u / \\partial U_x$ [px/px]'
	titles['dudUy'] = 'Displacement Derivative $\\partial u / \\partial U_y$ [px/px]'
	titles['dudVx'] = 'Displacement Derivative $\\partial u / \\partial V_x$ [px/px]'
	titles['dudVy'] = 'Displacement Derivative $\\partial u / \\partial V_y$ [px/px]'
	titles['v']     = 'Displacement $v$ [px]'
	titles['dvdU']  = 'Displacement Derivative $\\partial v / \\partial U$ [px/px]'
	titles['dvdV']  = 'Displacement Derivative $\\partial v / \\partial V$ [px/px]'
	titles['dvdUx'] = 'Displacement Derivative $\\partial v / \\partial U_x$ [px/px]'
	titles['dvdUy'] = 'Displacement Derivative $\\partial v / \\partial U_y$ [px/px]'
	titles['dvdVx'] = 'Displacement Derivative $\\partial v / \\partial V_x$ [px/px]'
	titles['dvdVy'] = 'Displacement Derivative $\\partial v / \\partial V_y$ [px/px]'
	D = {k: pl.empty( (Nf,Np) ) for k in keys}
	
	for fk in range(Nf):
		fn = fns[fk]
		fh = netcdf_file(fn)
		for key in keys:
			D[key][fk,:] = fh.variables[key][:,0,0]
	
	for key in keys:
		print(key)
		for pk in range(Np):
			name  = 'stats-%s-%d'%(key,pk)
			title = 'Pass %s: %s'%(pk,titles[key])
			plotHistogram(D[key][:,pk],name,title)

#fns = [ 'results/monteCarlo/pair-%d.nc'%(k+1) for k in range(Np) ]
#doPair(fns)

fns = [ 'results/monteCarlo/vectors-%d.nc'%(k+1) for k in range(Np) ]
doVector(fns)
