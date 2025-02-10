#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 2G
#SBATCH -c 1

#This script should be run from the results/cleanup_ASV_wise folder
#The file format of the output file differs from

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_sed_extr.txt
    for fname in autumn/cont_list_sed_extr*.txt
    do
        tail -n +2 $fname >> autumn/autumn_sed_extr.txt
    done

head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_sed_extr.txt
    for fname in spring/cont_list_sed_extr*.txt
    do
        tail -n +2 $fname >> spring/spring_sed_extr.txt
    done

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_sed_ntc.txt
    for fname in autumn/cont_list_sed_ntc*.txt
    do
        tail -n +2 $fname >> autumn/autumn_sed_ntc.txt
    done

head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_sed_ntc.txt
    for fname in spring/cont_list_sed_ntc*.txt
    do
        tail -n +2 $fname >> spring/spring_sed_ntc.txt
    done

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_wat_clst.txt
    for fname in autumn/cont_list_wat_clst*.txt
    do
        tail -n +2 $fname >> autumn/autumn_wat_clst.txt
    done

head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_wat_clst.txt
    for fname in spring/cont_list_wat_clst*.txt
    do
        tail -n +2 $fname >> spring/spring_wat_clst.txt
    done

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_wat_extr.txt
    for fname in autumn/cont_list_wat_extr*.txt
    do
        tail -n +2 $fname >> autumn/autumn_wat_extr.txt
    done

head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_wat_extr.txt
    for fname in spring/cont_list_wat_extr*.txt
    do
        tail -n +2 $fname >> spring/spring_wat_extr.txt
    done

head -n1 autumn/cont_list_sed_extr_autumn_2CVS4_.txt > autumn/autumn_wat_ntc.txt
    for fname in autumn/cont_list_wat_ntc*.txt
    do
        tail -n +2 $fname >> autumn/autumn_wat_ntc.txt
    done

head -n1 spring/cont_list_sed_extr_spring_S6_.txt > spring/spring_wat_ntc.txt
    for fname in spring/cont_list_wat_ntc*.txt
    do
        tail -n +2 $fname >> spring/spring_wat_ntc.txt
    done