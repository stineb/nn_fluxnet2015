#!/bin/bash

#for idx in `cat ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/test.txt`
for idx in `cat ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/sitelist_doprofile_fluxnet2015.txt`
do 
    cp  ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/submit_fVAR_nn_XXXXXX.sh  ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/submit_fVAR_nn_${idx}.sh
    sed -i -- s/XXXXX/${idx}/g  ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/submit_fVAR_nn_${idx}.sh
    echo "submitting fVAR script for site ${idx}"
    qsub  ${HOME}/sofun/utils_sofun/analysis_sofun/fluxnet2015/submit_fVAR_nn_${idx}.sh	
done
