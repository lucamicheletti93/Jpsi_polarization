#!/bin/bash

file="/media/luca/488AE2208AE20A70/PbPb_2015_Trees/run_list.txt"
while IFS= read line
do
	tree="/media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_$line.root"
	if [ -f "$tree" ]; then
	root -b -q create_trigger_response_function_from_data.C++\($line\)
	else
	echo "the file Tree_$line.root does not exist"
	fi
done <"$file"
