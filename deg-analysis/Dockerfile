FROM staphb/samtools:1.9

# ENV DST=/opt/samtools
# ENV ULOCAL=/usr/local/bin
ENV OUTPUT_DIR=output
ENV HTSEQ_FOLDER=htseq
ENV DEG_FOLDER=edgeR

USER root
RUN mkdir -p $HTSEQ_FOLDER
RUN mkdir -p $DEG_FOLDER

RUN apt-get clean -y && apt-get update -y && apt-get install -y rsync build-essential python-numpy python-matplotlib python-pysam python-htseq perl r-base
RUN mkdir -p $HTSEQ_FOLDER/$OUTPUT_DIR
RUN mkdir -p $DEG_FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/6-run-deg-analysis.sh"]
WORKDIR $HTSEQ_FOLDER/$OUTPUT_DIR

# CMD tail -f /dev/null