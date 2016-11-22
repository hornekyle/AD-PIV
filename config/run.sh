#!/bin/bash

base="/home/kyle/Work/Research/Branches/PIV/AD-PIV"
build="${base}/build"
template="${base}/config/template.cfg"
config="input.cfg"

export Npair=16

radii=(0.5 1.5 2.5)
disps=(5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0)
shears=(0.0 0.1 0.2 0.3 0.4 0.5)

for r in "${radii[@]}";
do
	for d in "${disps[@]}";
	do
		echo "Disp ${r}-${d}"
		export radius=${r}
		export disp=${d}
		export shear=0.0
		
		results="${base}/results/disp-${radius}"
		run="${disp}"
		mkdir -p "${results}/${run}"
		rm -f "${results}/${run}/*"
		cd "${results}/${run}"
		cat ${template} | envsubst > ${config}
		mpirun -np $(nproc) ${build}/process ${config} > output.log
		tar cfz "../d-${run}-data.tar.gz" "."
		cd ..
		rm -rf ${run}
	done
	
	for s in "${shears[@]}";
	do
		echo "Shear ${r}-${s}"
		export radius=${r}
		export disp=5.0
		export shear=${s}
		
		results="${base}/results/shear-${radius}"
		run="${shear}"
		mkdir -p "${results}/${run}"
		rm -f "${results}/${run}/*"
		cd "${results}/${run}"
		cat ${template} | envsubst > ${config}
		mpirun -np $(nproc) ${build}/process ${config} > output.log
		tar cfz "../s-${run}-data.tar.gz" "."
		cd ..
		rm -rf ${run}
	done
done
