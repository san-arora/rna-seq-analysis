#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/*/
do
  parent=$(basename "$sample")
  fastqc -o "$parent" -t 4 --extract "$sample"/*.gz
done