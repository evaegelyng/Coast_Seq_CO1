#!/bin/bash

for i in $(seq 1 20000) 
	do 
  		gwf run tab_COSQ_${i}
 	done
	echo "targets 1 to 20000 submitted"
