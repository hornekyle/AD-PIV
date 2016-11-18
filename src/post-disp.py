#!/usr/bin/env python

import matplotlib as mpl
#mpl.rcParams['font.family'] = 'sans'
#mpl.rcParams['font.size'] = 11
#mpl.rcParams['font.serif'] = 'Calibri'
mpl.rcParams['image.cmap'] = 'viridis'

import pylab as pl
from scipy.io import netcdf_file
import tarfile
import scipy.io

fext = 'png'
fdir = ''

radii = pl.linspace(0.5,2.5,2+1)
disps = pl.linspace(5.0,6.0,10+1)
Nd = disps.size

varNames = ['u','dudU','dudV','dudUx','dudUy','dudVx','dudVy','v','dvdU','dvdV','dvdUx','dvdUy','dvdVx','dvdVy']
passNames = ['true','map','map-shift','lsq-shift']
Np = len(passNames)

titles = {}
titles['u']     = 'Displacement\n$u$ [px]'
titles['dudU']  = 'Displacement Derivative\n$\\partial u / \\partial U$ [px/px]'
titles['dudV']  = 'Displacement Derivative\n$\\partial u / \\partial V$ [px/px]'
titles['dudUx'] = 'Displacement Derivative\n$\\partial u / \\partial U_x$ [px/1]'
titles['dudUy'] = 'Displacement Derivative\n$\\partial u / \\partial U_y$ [px/1]'
titles['dudVx'] = 'Displacement Derivative\n$\\partial u / \\partial V_x$ [px/1]'
titles['dudVy'] = 'Displacement Derivative\n$\\partial u / \\partial V_y$ [px/1]'
titles['v']     = 'Displacement\n$v$ [px]'
titles['dvdU']  = 'Displacement Derivative\n$\\partial v / \\partial U$ [px/px]'
titles['dvdV']  = 'Displacement Derivative\n$\\partial v / \\partial V$ [px/px]'
titles['dvdUx'] = 'Displacement Derivative\n$\\partial v / \\partial U_x$ [px/1]'
titles['dvdUy'] = 'Displacement Derivative\n$\\partial v / \\partial U_y$ [px/1]'
titles['dvdVx'] = 'Displacement Derivative\n$\\partial v / \\partial V_x$ [px/1]'
titles['dvdVy'] = 'Displacement Derivative\n$\\partial v / \\partial V_y$ [px/1]'

def getData(base):
	def getVarData(fn,var):
		def getMembers(tfh):
			names = sorted(tfh.getnames())
			names.remove('.')
			names.remove('./input.cfg')
			names.remove('./output.log')
			
			return [tfh.getmember(n) for n in names]

		def getPassCount(fl):
			fh = netcdf_file(fl)
			S = fh.variables['u'].shape
			fh.close()
			return S[0]
		
		tfh = tarfile.open(fn)
		members = getMembers(tfh)
		Nf = len(members)
		Np = getPassCount(tfh.extractfile(members[0]))
		
		data = pl.empty( (Nf,Np) )
		for k in range(len(members)):
			fl = tfh.extractfile(members[k])
			fh = netcdf_file(fl)
			data[k,:] = fh.variables[var][:,0,0]
			fh.close()
		tfh.close()
		return data
	
	# Find MC sample count
	fn = '%s/d-%3.1f-data.tar.gz'%(base,disps[0])
	data = getVarData(fn,'u')
	Ns = data.shape[0]
	# Create empty dictionary
	Data = {}
	# Pre-allocate arrays for storage
	for var in varNames:
		Data[var] = pl.empty( (Ns,Nd,Np) )
	# Fill arrays
	for disp_k in range(Nd):
		fn = '%s/d-%3.1f-data.tar.gz'%(base,disps[disp_k])
		for vn in varNames:
			data = getVarData(fn,vn)
			Data[vn][:,disp_k,:] = data[:,:]
	
	# Return Data[var][sample,displacment,pass]
	return Data

def plotStats(var,Data):
	def unzoom(ax,d,r):
		if d=='x':
			xl = ax.get_xlim()
			xr = xl[1]-xl[0]
			ax.set_xlim(xl[0]-xr*r,xl[1]+xr*r)
		if d=='y':
			yl = ax.get_ylim()
			yr = yl[1]-yl[0]
			ax.set_ylim(yl[0]-yr*r,yl[1]+yr*r)
	
	linesFormats = ['','g-o','c-s','m-^']
	plotTypes = ['mean','std']
	plotTypeNames = {'mean':'Mean Error','std':'Deviation'}
	for plotType in plotTypes:
		fig = pl.figure(figsize=(5,4),tight_layout={'rect':(0,0,1,0.93)})
		ax = fig.add_subplot(1,1,1)
		for pass_k in range(1,Np):
			x = disps
			if plotType=='mean':
				y = Data[var][:,:,pass_k].mean(0)-Data[var][:,:,0].mean(0)
			elif plotType=='std':
				y = Data[var][:,:,pass_k].std(0)
			ax.plot(x,y,linesFormats[pass_k],label=passNames[pass_k],mew=0)
		ax.set_xlim(disps.min(),disps.max())
		ax.set_xlabel('Displacement $u$ [px]')
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		unzoom(ax,'y',0.1)
		ax.set_ylabel('%s %s'%(plotTypeNames[plotType],titles[var]))
		ax.legend(loc='lower left', bbox_to_anchor=(-0.25,1.01),numpoints=3,
			frameon=False,ncol=3,borderpad=0.1,handletextpad=0.2,columnspacing=1.5)
		fig.savefig('%s/%s-%s.%s'%(fdir,var,plotType,fext))
		pl.close(fig)

def plotHistogram(ax,I,name,offset=0.0,scale=1.0,title='',Nb=None):
	if Nb==None:
		Nb = int(1.0*pl.sqrt(I.size))
	chance,edges = pl.histogram(I,Nb,density=True)
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

def plotHist(var,Data):
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout={'w_pad':0},figsize=(5,4))
		ax = fig.add_subplot(1,1,1)
		for disp_k in range(Nd):
			I = Data[var][:,disp_k,pass_k]
			plotHistogram(ax,I,titles[var],disps[disp_k],0.05)
		ax.set_xlim(disps.min()-0.1,disps.max()+0.1)
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		ax.set_xlabel('Displacement $u$ [px]')
		ax.set_ylabel('%s'%titles[var])
		fig.savefig('%s/%s-%d-hist.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

for rk in range(radii.size):
	base = './results/disp-%.1f'%radii[rk]
	fn = '%s/combined.mat'%base
	fdir = 'figures/disp-%.1f'%radii[rk]
	try:
		Data = scipy.io.loadmat(fn)
	except:
		Data = getData(base)
		scipy.io.savemat(fn,Data)


	for var in varNames:
		plotStats(var,Data)
		plotHist(var,Data)
