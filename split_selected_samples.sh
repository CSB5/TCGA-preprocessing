#/bin/bash

awk 'BEGIN {FS="-"} {if($4 ~ /^1/) print $0}' $1 > selected_normals.dat
awk 'BEGIN {FS="-"} {if($4 ~ /^0/) print $0}' $1 > excluded_tumors.dat

