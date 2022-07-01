################# BASE IMAGE ######################
FROM biocontainers/biocontainers:v1.1.0_cv2

################## METADATA ######################
LABEL base_image="biocontainers:v1.1.0_cv2"
LABEL version="8"
LABEL software="fastqc"
LABEL software.version="0.11.9"
LABEL about.summary="A quality control tool for high throughput sequence data."
LABEL about.home="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
LABEL about.documentation="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/"
LABEL about.license_file="https://www.gnu.org/copyleft/gpl.html"
LABEL about.license="SPDX:GPL-3.0"
LABEL extra.identifiers.biotools="fastqc"
LABEL about.tags="General"
LABEL extra.binaries="fastqc"


################## MAINTAINER ######################
LABEL maintainer="Alvin Chen <ychen@aperiomics.com>"

################## INSTALLATION ######################
USER root

ENV ZIP=fastqc_v0.11.9.zip
ENV URL=https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV FOLDER=FastQC
ENV DST=/opt/fastqc
ENV ULOCAL=/usr/local/bin
ENV OUTPUT_DIR=output

RUN mkdir -p $FOLDER
RUN mkdir -p $DST

RUN apt-get update && apt-get install -y openjdk-8-jre-headless && apt-get install -y rsync \
    wget --no-check-certificate $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER && chmod +x $DST/fastqc && \
    ln -s $DST/fastqc $ULOCAL/fastqc

COPY . /data

RUN mkdir $FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/1-run-fastqc.sh"]
WORKDIR $FOLDER/$OUTPUT_DIR

# CMD ["fastqc", "-o", "./",  "-t", "4", "--extract", "/data/dataset/1-Cas9-GFP/1-Ca9-GFP_S37_R1_001.fastq.gz"]
# CMD tail -f /dev/null
