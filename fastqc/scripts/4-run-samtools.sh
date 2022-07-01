#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/**/*
do
  parent=$(basename $(dirname "$sample"))
  fileName=$(basename "$sample")
  fileNameWithoutExt=$(basename "$sample" | cut -d. -f1)
  if [ ! -f "$parent"/"$fileNameWithoutExt".bam ]; then
    samtools view -S -h -b /data/hisat2/output/"$parent"/"$fileNameWithoutExt".sam > "$parent"/"$fileNameWithoutExt".bam
  fi
done