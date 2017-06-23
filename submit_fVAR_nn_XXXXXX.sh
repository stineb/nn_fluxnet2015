#!/bin/sh

## This tells the batch manager to limit the walltime for the job to XX hours, YY minutes and ZZ second
#PBS -l walltime=00:20:00

## This tells the batch manager to use XX nodes with YY cpus (total XX*YY cpus) and ZZ mb of memory per node.
#PBS -l select=1:ncpus=1:mem=1gb

## and use PP gb of memory.

module load R
module unload liblzma

## This jobs requires the Intel math kernel so we must load it at run time.

R CMD BATCH --no-save --no-restore '--args sitename="XXXXX" use_weights=FALSE use_fapar=FALSE nam_target="lue_obs_evi"' $HOME/sofun/utils_sofun/analysis_sofun/fluxnet2015/execute_nn_fVAR_fluxnet2015.R $HOME/sofun/utils_sofun/analysis_sofun/fluxnet2015/fvar_XXXXX.Rout

## This tells the batch manager to execute the program lazy from the examples
## directory of the users home directory.
