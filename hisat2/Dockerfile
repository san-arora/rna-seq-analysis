FROM biocontainers/biocontainers:vdebian-buster-backports_cv1
LABEL    software="hisat2" \ 
    base_image="biocontainers/biocontainers:vdebian-buster-backports_cv1" \ 
    container="hisat2" \ 
    about.summary="graph-based alignment of short nucleotide reads to many genomes" \ 
    about.home="https://ccb.jhu.edu/software/hisat2/index.shtml" \ 
    software.version="2.1.0-2-deb" \ 
    upstream.version="2.1.0" \ 
    version="1" \ 
    extra.identifiers.biotools="hisat2" \ 
    about.copyright=" 2011-2012, Ben Langmead <langmea@cs.jhu.edu>" \ 
    about.license="GPL-3+" \ 
    about.license_file="/usr/share/doc/hisat2/copyright" \ 
    extra.binaries="/usr/bin/hisat2,/usr/bin/hisat2-align-l,/usr/bin/hisat2-align-s,/usr/bin/hisat2-build,/usr/bin/hisat2-build-l,/usr/bin/hisat2-build-s,/usr/bin/hisat2-inspect,/usr/bin/hisat2-inspect-l,/usr/bin/hisat2-inspect-s" \ 
    about.tags=""

ENV DST=/opt/hisat2
ENV ULOCAL=/usr/local/bin
ENV OUTPUT_DIR=output
ENV FOLDER=hisat2

USER root
RUN mkdir -p $DST
RUN mkdir -p $FOLDER
COPY . /data/$FOLDER/

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get clean -y && apt-get update -y && apt-get install -y unzip && apt-get install -y python3 && apt-get install -y rsync
RUN unzip /data/$FOLDER/binary/hisat2-2.2.1.zip -d $DST
RUN ln -s $DST/hisat2-2.2.1/hisat2 $ULOCAL/hisat2
RUN ln -s $DST/hisat2-2.2.1/hisat2-align-l $ULOCAL/hisat2-align-l
RUN ln -s $DST/hisat2-2.2.1/hisat2-align-s $ULOCAL/hisat2-align-s
RUN ln -s $DST/hisat2-2.2.1/hisat2-build $ULOCAL/hisat2-build
RUN ln -s $DST/hisat2-2.2.1/hisat2-build-l $ULOCAL/hisat2-build-l
RUN ln -s $DST/hisat2-2.2.1/hisat2-build-s $ULOCAL/hisat2-build-s
RUN ln -s $DST/hisat2-2.2.1/hisat2-inspect $ULOCAL/hisat2-inspect
RUN ln -s $DST/hisat2-2.2.1/hisat2-inspect-l $ULOCAL/hisat2-inspect-l
RUN ln -s $DST/hisat2-2.2.1/hisat2-inspect-s $ULOCAL/hisat2-inspect-s

RUN mkdir $FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/3-run-hisat2.sh"]
WORKDIR $FOLDER/$OUTPUT_DIR

# CMD ["hisat2", "-x", "/data/reference-genome/grch37_snp_tran/genome_snp_tran", "-U", "/data/cutadapt/output/1-Ca9-GFP_S37_R1_001_cutadapt.fastq.gz", "-S", "1-Ca9-GFP_S37_R1_001_cutadapt.sam", "-p", "8", "--qc-filter", "--summary-file", "1-Ca9-GFP.hisat2.summary.txt"]
# CMD tail -f /dev/null
