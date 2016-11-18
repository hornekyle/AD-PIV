#!/bin/bash

base="/home/kyle/Work/Research/Branches/PIV/AD-PIV"
build="${base}/build"
template="${base}/config/template-disp.cfg"
config="input.cfg"

export Npair=64
export shear=0.0

radii=(0.5 1.5 2.5)
disps=(5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0)

for r in "${radii[@]}";
do
	for d in "${disps[@]}";
	do
		echo "${r}-${d}"
		export radius=${r}
		export disp=${d}
		
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
done
