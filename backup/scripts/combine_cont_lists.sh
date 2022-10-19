#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH -c 1

#This script should be run from the results/cleanup_ASV_wise folder

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_cont_list.txt
    for fname in autumn/*.txt
    do
        tail -n +2 $fname >> autumn/autumn_cont_list.txt
    done

#Something goes wrong in the spring list, so the file keeps growing

#head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_cont_list.txt
#    for fname in spring/*.txt
#    do
#        tail -n +2 $fname >> spring/spring_cont_list.txt
#    done