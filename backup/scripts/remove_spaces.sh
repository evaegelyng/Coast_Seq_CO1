#!/bin/bash

input_file=$1
output_file=$2

sed 's/ *$//g' $input_file > $output_file