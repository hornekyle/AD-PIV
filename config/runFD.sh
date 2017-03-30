#!/bin/bash

base="/home/kyle/Work/Research/Branches/PIV/AD-PIV"
build="${base}/build"
template="${base}/config/templateFD.cfg"
config="input.cfg"

export Npair=1
export radius=1.5

disps=(4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6)

for d in "${disps[@]}";
do
	echo "Disp ${d}"
	export disp=${d}
	
	results="${base}/results/finiteDifference"
	run="${disp}"
	mkdir -p "${results}/${run}"
	rm -f "${results}/${run}/*"
	cd "${results}/${run}"
	cat ${template} | envsubst > ${config}
	${build}/process ${config} > output.log
done
