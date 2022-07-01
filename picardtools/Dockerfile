FROM broadinstitute/picard
LABEL    software="picard-tools" \ 
    base_image="broadinstitute/picard" \ 
    container="picard-tools" \ 
    about.summary="Command line tools to manipulate SAM and BAM files" \ 
    about.home="http://broadinstitute.github.io/picard/" \ 
    software.version="2.18.25dfsg-2-deb" \ 
    upstream.version="2.18.25" \ 
    version="1" \ 
    about.copyright=" 2009-2018 The Broad Institute" \ 
    about.license="Expat" \ 
    about.license_file="/usr/share/doc/picard-tools/copyright" \ 
    extra.binaries="/usr/bin/PicardCommandLine,/usr/bin/picard-tools" \ 
    about.tags=""

ENV OUTPUT_DIR=output
ENV FOLDER=picardtools

USER root
RUN mkdir -p /data/$FOLDER
RUN apt-get update -y && apt-get install -y rsync
RUN mkdir /data/$FOLDER/$OUTPUT_DIR
ENTRYPOINT ["/bin/sh", "/data/scripts/5-run-picardtools.sh"]
WORKDIR /data/$FOLDER/$OUTPUT_DIR

# CMD ["java", "-jar", "/usr/picard/picard.jar", "SortSam", "I=/data/samtools/output/1-Ca9-GFP_S37_R1_001_bam.bam", "O=/data/picardtools/output/1-Ca9-GFP_sorted.bam", "SORT_ORDER=coordinate", "VALIDATION_STRINGENCY=SILENT", "CREATE_INDEX=true"]
# CMD tail -f /dev/null