#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/**/*
do
  parent=$(basename $(dirname "$sample"))
  fileName=$(basename "$sample")
  fileNameWithoutExt=$(basename "$sample" | cut -d. -f1)
  if [ ! -f "$parent"/"$fileNameWithoutExt"_sorted.bam ]; then
    java -jar /usr/picard/picard.jar SortSam I=/data/samtools/output/"$parent"/"$fileNameWithoutExt".bam O="$parent"/"$fileNameWithoutExt"_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
  fi
done