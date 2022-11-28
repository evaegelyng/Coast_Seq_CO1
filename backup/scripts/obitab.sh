#!/bin/bash
#SBATCH --partition normal 
#SBATCH --mem 196G         
#SBATCH -c 1               
#SBATCH -t 36:00:00        

#Run from main folder and activate mjolnir environment

obitab -o tmp/COSQ.new.fasta > tmp/COSQ_new.3.tab

