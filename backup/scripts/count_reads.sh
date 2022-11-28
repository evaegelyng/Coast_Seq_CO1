#!/bin/bash
#SBATCH --partition normal  
#SBATCH --mem 4G            
#SBATCH -c 1                
#SBATCH -t 04:00:00         

#Run from tmp folder

for d in W*; do
  [[ -d "$d" ]] && cd "$d" || continue
  f=".joined.fasta" 
  printf "$d\t$(cat ${f} | grep -c "^>")\n" >> ../filtered_read_counts.tsv
  cd -
done

for d in S*; do
  [[ -d "$d" ]] && cd "$d" || continue
  f=".joined.fasta" 
  printf "$d\t$(cat ${f} | grep -c "^>")\n" >> ../filtered_read_counts.tsv
  cd -
done