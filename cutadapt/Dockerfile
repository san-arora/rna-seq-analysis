FROM biocontainers/biocontainers:vdebian-buster-backports_cv1
LABEL    software="cutadapt" \ 
    base_image="biocontainers/biocontainers:vdebian-buster-backports_cv1" \ 
    container="cutadapt" \ 
    about.summary="Clean biological sequences from high-throughput sequencing reads (Python 3)" \ 
    about.home="http://pypi.python.org/pypi/cutadapt" \ 
    software.version="1.18-1-deb-py3" \ 
    upstream.version="1.18" \ 
    version="1" \ 
    about.copyright="2010-2015 Marcel Martin <marcel.martin@scilifelab.se>" \ 
    about.license="Expat" \ 
    about.license_file="/usr/share/doc/cutadapt/copyright" \ 
    about.tags=""

USER root
ENV DEBIAN_FRONTEND noninteractive
ENV OUTPUT_DIR=output
ENV ULOCAL=/usr/local/bin
ENV FOLDER=cutadapt

RUN mkdir -p $FOLDER

RUN apt-get update && apt-get install -y python3 python3-pip && apt-get install -y rsync
RUN python3 -m pip install --user --upgrade cutadapt
RUN ln -s $HOME/.local/bin/cutadapt $ULOCAL/cutadapt

COPY . /data/$FOLDER
RUN mkdir $FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/2-run-cutadapt.sh"]
WORKDIR $FOLDER/$OUTPUT_DIR

# CMD ["cutadapt", "-a", "AAAAAAAAAA",  "-q", "10", "-m", "20", "-u", "12", "-o", "1-Ca9-GFP_S37_R1_001_cutadapt.fastq.gz", "/data/dataset/1-Cas9-GFP/1-Ca9-GFP_S37_R1_001.fastq.gz"]
# CMD tail -f /dev/null