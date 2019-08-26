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

radii = pl.linspace(0.5,2.5,2+1)
disps = pl.linspace(5.0,6.0,10+1)

Nd = disps.size

SP = 50

#varNames = ['u','dudU','dudV','dudUx','dudUy','dudVx','dudVy','v','dvdU','dvdV','dvdUx','dvdUy','dvdVx','dvdVy']
varNames = ['u','dudU','dudV','v','dvdU','dvdV']
passNames = ['true','map','map-shift','lsq-shift']
Np = len(passNames)

titles = {}
titles['u']     = 'Displacement $u$ $/$ px'
titles['dudU']  = 'Displacement Derivative $\\partial_U (u)$ $/$ 1'
titles['dudV']  = 'Displacement Derivative $\\partial_V (u)$ $/$ 1'
#titles['dudUx'] = 'Displacement Derivative $\\partial_{U_x} (u)$ $/$ px'
#titles['dudUy'] = 'Displacement Derivative $\\partial_{U_y} (u)$ $/$ px'
#titles['dudVx'] = 'Displacement Derivative $\\partial_{V_x} (u)$ $/$ px'
#titles['dudVy'] = 'Displacement Derivative $\\partial_{V_y} (u)$ $/$ px'
titles['v']     = 'Displacement $v$ [px]'
titles['dvdU']  = 'Displacement Derivative $\\partial_U (v)$ $/$ 1'
titles['dvdV']  = 'Displacement Derivative $\\partial_V (v)$ $/$ 1'
#titles['dvdUx'] = 'Displacement Derivative $\\partial_{U_x} (v)$ $/$ px'
#titles['dvdUy'] = 'Displacement Derivative $\\partial_{U_y} (v)$ $/$ px'
#titles['dvdVx'] = 'Displacement Derivative $\\partial_{V_x} (v)$ $/$ px'
#titles['dvdVy'] = 'Displacement Derivative $\\partial_{V_y} (v)$ $/$ px'

shortTitles = {}
shortTitles['u']     = 'u'
shortTitles['dudU']  = '\\partial_U (u)'
shortTitles['dudV']  = '\\partial_V (u)'
#shortTitles['dudUx'] = '\\partial_{U_x} (u)'
#shortTitles['dudUy'] = '\\partial_{U_y} (u)'
#shortTitles['dudVx'] = '\\partial_{V_x} (u)'
#shortTitles['dudVy'] = '\\partial_{V_y} (u)'
shortTitles['v']     = 'v'
shortTitles['dvdU']  = '\\partial_U (v)'
shortTitles['dvdV']  = '\\partial_V (v)'
#shortTitles['dvdUx'] = '\\partial_{U_x} (v)'
#shortTitles['dvdUy'] = '\\partial_{U_y} (v)'
#shortTitles['dvdVx'] = '\\partial_{V_x} (v)'
#shortTitles['dvdVy'] = '\\partial_{V_y} (v)'

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

def plotStats(var,Data):
	linesFormats = ['','C0-o','C1-s','C2-^']
	plotTypes = ['mean','std']
	plotTypeNames = {'mean':'Mean Error','std':'Deviation'}
	for plotType in plotTypes:
		fig = pl.figure(figsize=(4,3),tight_layout={'rect':(0,0,1,0.93)})
		ax = fig.add_subplot(1,1,1)
		for pass_k in range(1,Np):
			x = disps
			if plotType=='mean':
				y = Data[var][:,:,pass_k].mean(0)-Data[var][:,:,0].mean(0)
			elif plotType=='std':
				y = Data[var][:,:,pass_k].std(0)
			ax.plot(x,y,linesFormats[pass_k],label=passNames[pass_k],mew=0)
		ax.set_xlim(disps.min(),disps.max())
		ax.set_xlabel('Correct Displacement $U$ [px]')
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		unzoom(ax,'y',0.1)
		ax.set_ylabel('%s\n%s'%(plotTypeNames[plotType],titles[var]))
		ax.legend(loc='lower left', bbox_to_anchor=(-0.15,1.01),numpoints=3,
			frameon=False,ncol=3,borderpad=0.1,handletextpad=0.2,columnspacing=1.5)
		fig.savefig('%s/%s-%s.%s'%(fdir,var,plotType,fext))
		pl.close(fig)

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

def plotHist(var,Data):
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout=True,figsize=(4,3))
		ax = fig.add_subplot(1,1,1)
		for disp_k in range(Nd):
			I = Data[var][:,disp_k,pass_k]
			plotHistogram(ax,I,titles[var],disps[disp_k],0.08)
			ax.errorbar(disps[disp_k],I.mean(),I.std(),fmt='k.')
		ax.set_xlim(disps.min()-0.1,disps.max()+0.1)
		ax.set_xticks(pl.linspace(5.0,6.0,5+1))
		ax.set_xlabel('Correct Displacement $U$ $/$ px')
		ax.set_ylim( computeBounds(Data[var][:,:,pass_k]) )
		unzoom(ax,'y',0.2)
		ax.set_ylabel('%s'%titles[var])
		fig.savefig('%s/%s-%d-hist.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

def plotViolin(var,Data):
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout=True,figsize=(4,3))
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
		ax.set_xlabel('Correct Displacement $U$ $/$ px')
		ax.set_ylim( computeBounds(Data[var][:,:,pass_k]) )
		unzoom(ax,'y',0.2)
		ax.set_ylabel('%s'%titles[var])
		fig.savefig('%s/%s-%d-violin.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

def plotConvergence(var,Data):
	Ns = Data[var][:,0,0].size
	for pass_k in range(1,Np):
		fig = pl.figure(tight_layout={'h_pad':0,'rect':(0,0,1,0.95)},figsize=(4,3))
		ax1 = fig.add_subplot(2,1,1)
		ax2 = fig.add_subplot(2,1,2)
		i = pl.linspace(1,Ns+1,Ns)
		m = pl.zeros( (Nd,Ns) )
		s = pl.zeros( (Nd,Ns) )
		I = Data[var][:,:,pass_k]
		for k in range(2,Ns):
			m[:,k] = I[:k,:].mean(0)
			s[:,k] = I[:k,:].std(0)
		for disp_k in range(Nd):
			m[disp_k,:] = m[disp_k,:]-m[disp_k,-1]
			s[disp_k,:] = s[disp_k,:]-s[disp_k,-1]
		m0 = 3.0*abs(m).mean()
		s0 = 5.0*abs(s).mean()
		for disp_k in range(Nd):
			ax1.plot(i[2::SP],m[disp_k,2::SP]/m0,'C0-')
			ax2.plot(i[2::SP],s[disp_k,2::SP]/s0,'C1-')
		
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
		ax2.set_xlabel(r'Monte Carlo Iterations $k$ $/$ \#')
		#fig.suptitle('Normalized %s'%titles[var])
		ax1.set_ylabel('Mean\n$n( \\mu_{%s} )$ $/$ 1'%shortTitles[var])
		ax2.set_ylabel('Deviation\n$n( \\delta_{%s} )$ $/$ 1'%shortTitles[var])
		fig.savefig('%s/%s-%d-conv.%s'%(fdir,var,pass_k,fext))
		pl.close(fig)

def plotCorrelation(Data):
	for var in ['u','v']:
		for dim in ['U','V']:
			for kp in range(1,Np):
				Ns = Data['d%sd%s'%(var,dim)][:,0,kp].shape[0]
				
				e = pl.empty( (Nd,Ns) )
				d = pl.empty( (Nd,Ns) )
				s = pl.empty( (Nd,Ns) )
				order = pl.array(range(d.size))
				shuffle(order)
				
				fig = pl.figure(figsize=(4,3),tight_layout=True)
				ax  = fig.add_subplot(1,1,1)
				for kd in range(Nd):
					t = 0.0
					if var=='u': t = disps[kd]
					e[kd,:] = Data[var][:,kd,kp]-t
					d[kd,:] = Data['d%sd%s'%(var,dim)][:,kd,kp]
					s[kd,:] = pl.ones(Ns)*disps[kd]
				e = e.flatten()[order]
				d = d.flatten()[order]
				s = s.flatten()[order]
				
				sp = ax.scatter(d[::SP],e[::SP],c=s[::SP],marker='+')
				ax.set_xlabel(titles['d%sd%s'%(var,dim)])
				ax.set_ylabel('Displacement Error $\\epsilon_%s$ $/$ px'%var)
				fig.colorbar(sp,ax=ax,label='Displacement $u$ $/$ px')
				fig.savefig('%s/d%sd%s-%d-corr.%s'%(fdir,var,dim,kp,fext),dpi=300)
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
		plotStats(var,Data)
		plotHist(var,Data)
		plotViolin(var,Data)
		plotConvergence(var,Data)
	plotCorrelation(Data)
