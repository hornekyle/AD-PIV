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
import tarfile
import scipy.io
from random import shuffle

fext = 'svg'
fdir = ''

#radii = pl.linspace(0.5,2.5,2+1)
radii = pl.array([0.5])
disps = pl.linspace(5.0,6.0,10+1)

Nd = disps.size

SP = 50

#varNames = ['u','dudU','dudV','dudUx','dudUy','dudVx','dudVy','v','dvdU','dvdV','dvdUx','dvdUy','dvdVx','dvdVy']
varNames = ['u','dudU','dudV','v','dvdU','dvdV']
passNames = ['true','map','map-shift','lsq-shift']
Np = len(passNames)

titles = {}
titles['u']     = 'Displacement $u$ [px]'
titles['dudU']  = 'Displacement Derivative $\\partial u / \\partial U$ [px/px]'
titles['dudV']  = 'Displacement Derivative $\\partial u / \\partial V$ [px/px]'
#titles['dudUx'] = 'Displacement Derivative $\\partial u / \\partial U_x$ [px/1]'
#titles['dudUy'] = 'Displacement Derivative $\\partial u / \\partial U_y$ [px/1]'
#titles['dudVx'] = 'Displacement Derivative $\\partial u / \\partial V_x$ [px/1]'
#titles['dudVy'] = 'Displacement Derivative $\\partial u / \\partial V_y$ [px/1]'
titles['v']     = 'Displacement $v$ [px]'
titles['dvdU']  = 'Displacement Derivative $\\partial v / \\partial U$ [px/px]'
titles['dvdV']  = 'Displacement Derivative $\\partial v / \\partial V$ [px/px]'
#titles['dvdUx'] = 'Displacement Derivative $\\partial v / \\partial U_x$ [px/1]'
#titles['dvdUy'] = 'Displacement Derivative $\\partial v / \\partial U_y$ [px/1]'
#titles['dvdVx'] = 'Displacement Derivative $\\partial v / \\partial V_x$ [px/1]'
#titles['dvdVy'] = 'Displacement Derivative $\\partial v / \\partial V_y$ [px/1]'

shortTitles = {}
shortTitles['u']     = 'u'
shortTitles['dudU']  = '\\partial u / \\partial U'
shortTitles['dudV']  = '\\partial u / \\partial V'
#shortTitles['dudUx'] = '\\partial u / \\partial U_x'
#shortTitles['dudUy'] = '\\partial u / \\partial U_y'
#shortTitles['dudVx'] = '\\partial u / \\partial V_x'
#shortTitles['dudVy'] = '\\partial u / \\partial V_y'
shortTitles['v']     = 'v'
shortTitles['dvdU']  = '\\partial v / \\partial U'
shortTitles['dvdV']  = '\\partial v / \\partial V'
#shortTitles['dvdUx'] = '\\partial v / \\partial U_x'
#shortTitles['dvdUy'] = '\\partial v / \\partial U_y'
#shortTitles['dvdVx'] = '\\partial v / \\partial V_x'
#shortTitles['dvdVy'] = '\\partial v / \\partial V_y'

def unzoom(ax,d,r):
		if d=='x':
			xl = ax.get_xlim()
			xr = xl[1]-xl[0]
			ax.set_xlim(xl[0]-xr*r,xl[1]+xr*r)
		if d=='y':
			yl = ax.get_ylim()
			yr = yl[1]-yl[0]
			ax.set_ylim(yl[0]-yr*r,yl[1]+yr*r)

# Return Data[var][sample,displacment,pass]
def getData(base):
	def getVarData(fn,var):
		fh = netcdf_file(fn)
		Np = fh.variables[var].shape[1]
		Nf = fh.variables[var].shape[0]
		data = pl.empty( (Nf,Np) )
		data[:,:] = fh.variables[var][:,:,0,0]
		fh.close()
		return data
	
	# Find MC sample count
	fn = '%s/%3.1f/vectors.nc'%(base,disps[0])
	data = getVarData(fn,'u')
	Ns = data.shape[0]
	# Create empty dictionary
	Data = {}
	# Pre-allocate arrays for storage
	for var in varNames:
		Data[var] = pl.empty( (Ns,Nd,Np) )
	# Fill arrays
	for disp_k in range(Nd):
		fn = '%s/%3.1f/vectors.nc'%(base,disps[disp_k])
		for vn in varNames:
			data = getVarData(fn,vn)
			Data[vn][:,disp_k,:] = data[:,:]
	return Data

def histogram(I,Nb):
	L = pl.sort(I)
	dL = int(L.size/Nb)
	edges = pl.empty(Nb+1)
	chance = pl.empty(Nb)
	for k in range(Nb):
		edges[k] = L[k*dL]
	edges[-1] = L[-1]
	for k in range(Nb):
		chance[k] = dL/(edges[k+1]-edges[k])
	return chance,edges
	

def plotHistogram(ax,I,name,offset=0.0,scale=1.0,title='',Nb=None):
	#if Nb==None:
		#Nb = int(1.0*pl.sqrt(I.size))
	#chance,edges = pl.histogram(I,Nb,density=True)
	chance,edges = histogram(I,25)
	chance = scale*(chance/chance.max())
	edgesLeft  = edges[:-1]
	edgesRight = edges[1:]
	centers = (edgesLeft+edgesRight)/2.0
	heights = edgesRight-edgesLeft
	
	smap = mpl.cm.ScalarMappable()
	smap.set_array(chance)
	colors = smap.to_rgba(chance)
	widths = chance
	left = offset-chance/2
	ax.barh(centers,widths,heights,left,lw=0,align='center',color=colors)

def computeBounds(I):
	chance,edges = pl.histogram(I,int(1.0*pl.sqrt(I.size)),density=True)
	mx = (edges[1:]+edges[:-1])/2.0
	dx = edges[1:]-edges[:-1]
	y = chance
	xl = mx.min()
	xh = mx.max()
	for k in range(y.size):
		if (y[:k]*dx[:k]).sum()>=0.01:
			xl = mx[k]
			break
	for k in range(1,y.size):
		if (y[-k:]*dx[-k:]).sum()>=0.01:
			xh = mx[-k]
			break
	return (xl,xh)

def plotViolin(var,Data):
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout=True,figsize=(4,4))
		ax = fig.add_subplot(1,1,1)
		for disp_k in range(Nd):
			I = Data[var][:,disp_k,pass_k]
			ax.errorbar(disps[disp_k],I.mean(),I.std(),fmt='k.')
		violin_parts = ax.violinplot(Data[var][:,:,pass_k],disps,vert=True,widths=0.075,showextrema=False)
		for pc in violin_parts['bodies']:
			pc.set_alpha(0.5)
			pc.set_facecolor('C0')
			pc.set_edgecolor('')
		ax.set_xlim(disps.min()-0.1,disps.max()+0.1)
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		ax.set_xlabel('Displacement $u$ [px]')
		ax.set_ylim( computeBounds(Data[var][:,:,pass_k]) )
		unzoom(ax,'y',0.2)
		ax.set_ylabel('%s'%titles[var])
		fig.savefig('%s/%s-%d-violin.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

# Plot all data series
for rk in range(radii.size):
	base = './results/disp-%.1f'%radii[rk]
	fn = '%s/combined.mat'%base
	fdir = 'figures/disp-%.1f'%radii[rk]
	# Pre-cache combined data if not done
	try:
		Data = scipy.io.loadmat(fn)
	except:
		Data = getData(base)
		scipy.io.savemat(fn,Data)
	for var in varNames:
		print(rk,var)
		plotViolin(var,Data)
