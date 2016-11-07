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

lyxDoc = r"""
#LyX 2.2 created this file. For more info see http://www.lyx.org/
\lyxformat 474
\begin_document
\begin_header
\textclass article
\use_default_options true
\maintain_unincluded_children false
\language english
\language_package default
\inputencoding auto
\fontencoding global
\font_roman default
\font_sans default
\font_typewriter default
\font_math auto
\font_default_family default
\use_non_tex_fonts false
\font_sc false
\font_osf false
\font_sf_scale 100
\font_tt_scale 100
\graphics default
\default_output_format default
\output_sync 0
\bibtex_command default
\index_command default
\paperfontsize default
\use_hyperref false
\papersize default
\use_geometry false
\use_package amsmath 1
\use_package amssymb 1
\use_package cancel 1
\use_package esint 1
\use_package mathdots 1
\use_package mathtools 1
\use_package mhchem 1
\use_package stackrel 1
\use_package stmaryrd 1
\use_package undertilde 1
\cite_engine basic
\cite_engine_type default
\biblio_style plain
\use_bibtopic false
\use_indices false
\paperorientation portrait
\suppress_date false
\justification true
\use_refstyle 1
\index Index
\shortcut idx
\color #008000
\end_index
\secnumdepth 3
\tocdepth 3
\paragraph_separation indent
\paragraph_indentation default
\quotes_language english
\papercolumns 1
\papersides 1
\paperpagestyle default
\tracking_changes false
\output_changes false
\html_math_output 0
\html_css_as_file 0
\html_be_strict false
\end_header

\begin_body

\begin_layout Standard

%(figures)s

\end_layout

\end_body
\end_document
"""

lyxFigure = r"""
\begin_inset Float figure
wide false
sideways false
status open

\begin_layout Plain Layout
\align center
\begin_inset Graphics
	filename %(fileName)s

\end_inset

\end_layout

\begin_layout Plain Layout
\begin_inset Caption Standard

\begin_layout Plain Layout
\begin_inset ERT
status open

\begin_layout Plain Layout

%(caption)s
\end_layout

\end_inset

\begin_inset CommandInset label
LatexCommand label
name "fig:%(label)s"

\end_inset

\end_layout

\end_inset


\end_layout

\end_inset
"""

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

#==================#
#= Plot Image Set =#
#==================#

def doSet(setName):
	doPair('results/%s/pair-1.nc'%setName,'figures/%s-pair'%setName)
	doVector('results/%s/vector-1-[1,1|1].nc'%setName,'figures/%s-vector'%setName)
	doMap('results/%s/map-1-[1,1|1].nc'%setName,'figures/%s-map'%setName)

#================#
#= Run Commands =#
#================#

sets = [ 'set-%d'%k for k in range(Nsets) ]

for s in sets:
	doSet(s)
	
	figures = ''
	for k in capkeys:
		D = {}
		D['fileName'] = k
		D['caption'] = captions[k].replace('\\','\\backslash\n')
		D['label'] = k.split('/')[1].split('.')[0]
		
		figures = '%s\n\n%s'%(figures,lyxFigure%D)

	fh = open('%s.lyx'%s,'w')
	fh.write(lyxDoc%{'figures':figures})
	fh.close()
	
	capkeys = []
	captions = {}

pl.show()
