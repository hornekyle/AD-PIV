project: AD-PIV
project_dir: ./src
source: true
exclude: plplotlib.f90
display: public
graph: true
output_dir: ./doc
author: Dr. Kyle Horne, Rodolfo Grullon Varela
email: hornek@uwplatt.edu
summary: AD-PIV
display: public

Automatic differentiation Particle-image Velocimetry

Rodolfo needs co-author credit on both papers
Try to get Mike Black's help on Bayesian paper, and credit him as co-author

Image generation mode paper figures
	* Automatic differentiation (diagrams)
		- Expression decomposition
		- Chain rule propagation
	* Demonstration of AD applied to PIV (plots)
		- Image pair and derivatives
		- Correlation map/derivative fields with key derivatives
		- Vector pair and derivatives
		- Window shifting example at vector level
	* Monte Carlo results of single settings (plots)
		- Statistical convergence of MC method (look at GUM)
		- Average and standard deviation images
		- Histograms of u,dudU,dudV, etc.
	* Monte Carlo results for range of displacements (plots)
		- 4.5-5.5 px displacement (11 samples)
		- Plot mean indicated vs true displacement with error bars (maybe violin plots?)
		- Plot mean indicated vs true sensitivities with error bars (maybe violin plots?)
	* Monte Carlo results for range of shears (plots)
		- 0.0-0.5 px/px shear rates (6 samples)
		- 5.5 px displacement
		- Plot mean indicated vs true displacement with error bars (maybe violin plots?)
		- Plot mean indicated vs true sensitivities with error bars (maybe violin plots?)

Full differentiation mode todo
	* Look more closely at relationship between dudU, dudI, dIdU
	* Preliminary treatment of pixel uncertainty using Bayesian statistics
	* Actual uncertainty recommendations for a basic test case
	* Evaluation of the uncertainty effectiveness using Monte Carlo

Full differentiation mode figures
	* Automatic differentiation (diagrams)
		- 3D derivative space
		- Pixel seeding
		- Explanation of data loss at edges
	* Simple image pair (plots)
		- Comparison of dIdU and dudI
		- Comparison of dIdV and dvdI
	* Bayesian pattern statistics representation (plots)
	* Histograms of vectors, with comparison to predicted uncertainties (plots)
