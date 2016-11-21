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

#radii = pl.linspace(0.5,2.5,2+1)
radii = pl.array([1.5])
#disps = pl.linspace(5.0,6.0,10+1)
disps = pl.linspace(5.0,6.0,5+1)
#shears = pl.linspace(0.0,0.5,5+1)

Nd = disps.size

varNames = ['u','dudU','dudV','dudUx','dudUy','dudVx','dudVy','v','dvdU','dvdV','dvdUx','dvdUy','dvdVx','dvdVy']
passNames = ['true','map','map-shift','lsq-shift']
Np = len(passNames)

titles = {}
titles['u']     = 'Displacement $u$ [px]'
titles['dudU']  = 'Displacement Derivative $\\partial u / \\partial U$ [px/px]'
titles['dudV']  = 'Displacement Derivative $\\partial u / \\partial V$ [px/px]'
titles['dudUx'] = 'Displacement Derivative $\\partial u / \\partial U_x$ [px/1]'
titles['dudUy'] = 'Displacement Derivative $\\partial u / \\partial U_y$ [px/1]'
titles['dudVx'] = 'Displacement Derivative $\\partial u / \\partial V_x$ [px/1]'
titles['dudVy'] = 'Displacement Derivative $\\partial u / \\partial V_y$ [px/1]'
titles['v']     = 'Displacement $v$ [px]'
titles['dvdU']  = 'Displacement Derivative $\\partial v / \\partial U$ [px/px]'
titles['dvdV']  = 'Displacement Derivative $\\partial v / \\partial V$ [px/px]'
titles['dvdUx'] = 'Displacement Derivative $\\partial v / \\partial U_x$ [px/1]'
titles['dvdUy'] = 'Displacement Derivative $\\partial v / \\partial U_y$ [px/1]'
titles['dvdVx'] = 'Displacement Derivative $\\partial v / \\partial V_x$ [px/1]'
titles['dvdVy'] = 'Displacement Derivative $\\partial v / \\partial V_y$ [px/1]'

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
	
	return Data

def plotStats(var,Data):
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
		ax.set_ylabel('%s\n%s'%(plotTypeNames[plotType],titles[var]))
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
		fig = pl.figure(tight_layout=True,figsize=(5,4))
		ax = fig.add_subplot(1,1,1)
		for disp_k in range(Nd):
			I = Data[var][:,disp_k,pass_k]
			plotHistogram(ax,I,titles[var],disps[disp_k],0.05)
			ax.errorbar(disps[disp_k],I.mean(),I.std(),fmt='k.')
		unzoom(ax,'y',0.1)
		ax.set_xlim(disps.min()-0.1,disps.max()+0.1)
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		ax.set_xlabel('Displacement $u$ [px]')
		ax.set_ylabel('%s'%titles[var])
		fig.savefig('%s/%s-%d-hist.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

def plotConvergence(var,Data):
	Ns = Data[var][:,0,0].size
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout={'h_pad':0,'rect':(0,0,1,0.95)},figsize=(5,4))
		ax1 = fig.add_subplot(2,1,1)
		ax2 = fig.add_subplot(2,1,2)
		for disp_k in range(Nd):
			I = Data[var][:,disp_k,pass_k]
			i = pl.linspace(1,Ns+1,Ns)
			m = pl.empty(Ns)
			d = pl.empty(Ns)
			m[-1] = I.mean()
			d[-1] = I.std()
			for k in range(2,Ns-1):
				m[k] = I[:k].mean()-m[-1]
				d[k] = I[:k].std()-d[-1]
			m[-1] = 0.0
			d[-1] = 0.0
			m /= abs(m[2:]).max()
			d /= abs(d[2:]).max()
			ax1.plot(i[2:],m[2:],'g-')
			ax2.plot(i[2:],d[2:],'c-')
		
		ax1.axhline(0.0,color='k',linestyle='--')
		ax2.axhline(0.0,color='k',linestyle='--')
		ax1.set_xlim(1,Ns)
		ax2.set_xlim(1,Ns)
		ax1.set_xticks([])
		ax1.set_ylim(-1.1,1.1)
		ax2.set_ylim(-1.1,1.1)
		ax1.set_yticks([-1,0,1])
		ax1.set_yticklabels(['','',''])
		ax2.set_yticks([-1,0,1])
		ax2.set_yticklabels(['','',''])
		ax2.set_xlabel('Monte Carlo Iterations $k$ [#]')
		fig.suptitle('Normalized %s'%titles[var])
		ax1.set_ylabel('Mean [-]')
		ax2.set_ylabel('Deviation [-]')
		fig.savefig('%s/%s-%d-conv.%s'%(fdir,var,pass_k,fext))
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
		plotStats(var,Data)
		plotHist(var,Data)
		plotConvergence(var,Data)
