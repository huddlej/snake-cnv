#!/bin/bash

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 repeats.out > repeats.bed"
    exit 1
fi

sed '1,3d' $1 | awk 'OFS="\t" { if($9 == "C") { strand="-" } else { strand = $9 } print $5,$6,$7,$11,(100-$2)*10,strand,$10 }'
