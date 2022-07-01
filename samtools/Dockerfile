# FROM biocontainers/biocontainers:vdebian-buster-backports_cv1
FROM staphb/samtools:1.9
LABEL    software="samtools" \ 
    base_image="biocontainers/biocontainers:vdebian-buster-backports_cv1" \ 
    container="samtools" \ 
    about.summary="processing sequence alignments in SAM, BAM and CRAM formats" \ 
    about.home="http://www.htslib.org/" \ 
    software.version="1.9-4-deb" \ 
    upstream.version="1.9" \ 
    version="1" \ 
    extra.identifiers.biotools="samtools" \ 
    about.copyright=" 20082016, Genome Research Ltd. (GRL)" \ 
    about.license="Expat" \ 
    about.license_file="/usr/share/doc/samtools/copyright" \ 
    extra.binaries="/usr/bin/ace2sam,/usr/bin/blast2sam.pl,/usr/bin/bowtie2sam.pl,/usr/bin/export2sam.pl,/usr/bin/interpolate_sam.pl,/usr/bin/maq2sam-long,/usr/bin/maq2sam-short,/usr/bin/md5fa,/usr/bin/md5sum-lite,/usr/bin/novo2sam.pl,/usr/bin/plot-bamstats,/usr/bin/psl2sam.pl,/usr/bin/sam2vcf.pl,/usr/bin/samtools,/usr/bin/samtools.pl,/usr/bin/seq_cache_populate.pl,/usr/bin/soap2sam.pl,/usr/bin/varfilter.py,/usr/bin/wgsim,/usr/bin/wgsim_eval.pl,/usr/bin/zoom2sam.pl" \ 
    about.tags="field::biology, field::biology:bioinformatics, implemented-in::c,:commandline, network::client, protocol::ftp, protocol::http,:program, scope::utility, uitoolkit::ncurses, use::analysing,:calculating,:filtering, works-with::biological-sequence" 

# ENV DST=/opt/samtools
# ENV ULOCAL=/usr/local/bin
ENV OUTPUT_DIR=output
ENV FOLDER=samtools

USER root
RUN mkdir -p $FOLDER
# COPY . /data/$FOLDER/
# RUN mkdir -p $DST

# ENV DEBIAN_FRONTEND noninteractive
RUN apt-get clean -y && apt-get update -y && apt-get install -y rsync
# RUN apt-get clean -y && apt-get update -y && apt-get install -y bzip2 make gcc libz-dev libncurses-dev liblzma-dev
# RUN apt-get install -y autoconf automake perl zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev
# RUN apt-get clean -y && apt-get update -y && apt-get install -y make gcc libbz2-dev
# RUN make -C /data/samtools-install/samtools-1.8
RUN mkdir $FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/4-run-samtools.sh"]
WORKDIR $FOLDER/$OUTPUT_DIR

# CMD ["samtools", "view", "-S", "-h", "-b", "/data/hisat2/output/1-Ca9-GFP_S37_R1_001_cutadapt.sam", ">", "/data/samtools/output/1-Ca9-GFP_S37_R1_001_bam.bam"]
# CMD tail -f /dev/null