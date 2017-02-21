#!/bin/bash

#           A line beginning with #PBS is a PBS directive.
#           PBS directive must come first; any directives
# 	   after the firse executable statement are ignored.

#PBS -N zz

#Specify the number of nodes requested and the number 
#of processor per node

#PBS -l nodes=5:ppn=1
#
###PBS -o stdout_file
###PBS -e stderr_file
#               Specify the maximum cpu and wall clock time. The wall 
# 	   	clock time should take possible queue waiting time into account
#               Format:   hhhh:mm:ss hours:minutes:seconds
#PBS -l   cput=4:00:00
#PBS -l  walltime=8:00:00
#		Specify the maximun amount of physical memory required per process.
#                   kb for kilobytes, mb for megabytes, gb for gigabytes
###PBS -l pmem=512mb
#######################################
#		                      #
#  Output some useful job information #
#				      #
#######################################
#NCPU='wc -l < $PBS_NODEFILE'
#Please Set INPU env variable
#INPUT=
echo -------------------------------------------
echo 'This Job is allocated on '${NCPU}' cpu(s)'
echo 'Job is running on nodes(s): '
echo -------------------------------------------
#
cat $PBS_NODEFILE
WDIR=/scratch/$USER.$PBS_JOBID
cd $WDIR
pwd
mpirun_rsh -rsh -np $PBS_NP -hostfile $PBS_NODEFILE /home/stud6/src3/provaMatteo/jacobi
