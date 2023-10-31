#!/bin/bash

for file in $(ls ./)
do
	content=$(less $file)
	echo "$content"\n >> recap_out.txt
done

