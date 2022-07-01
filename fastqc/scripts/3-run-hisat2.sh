#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/**/*
do
  parent=$(basename $(dirname "$sample"))
  fileName=$(basename "$sample")
  fileNameWithoutExt=$(basename "$sample" | cut -d. -f1)
  if [ ! -f "$parent"/"$fileNameWithoutExt".sam ]; then
    hisat2 -x /data/reference-genome/grch37_snp_tran/genome_snp_tran -U /data/cutadapt/output/"$parent"/"$fileName" -S "$parent"/"$fileNameWithoutExt".sam -p 8 --qc-filter --summary-file "$parent"/"$fileNameWithoutExt".hisat2.summary.txt
  fi
done