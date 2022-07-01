#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/**/*
do
  parent=$(basename $(dirname "$sample"))
  fileName=$(basename "$sample")
  if [ ! -f "$parent"/"$fileName" ]; then
    cutadapt -a AAAAAAAAAA -q 10 -m 20 -u 12 -o "$parent"/"$fileName" "$sample"
  fi
done