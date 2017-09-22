#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: CRUK Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Mode: BY_SAMPLE
version="1.0.0"

# Directory structure required for pipeline
#
# /data
# └── results
#     └── seqId
#         ├── panel1
#         │   ├── sample1
#         │   ├── sample2
#         │   └── sample3
#         └── panel2
#             ├── sample1
#             ├── sample2
#             └── sample3
#
# Script 1 runs in sample folder, requires fastq files split by lane

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/CRUK/CRUK-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapters
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 35 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq.gz \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq.gz \
    "$read1Fastq" \
    "$read2Fastq"

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

done

#Print QC metrics
echo -e "RawSequenceQuality" > "$seqId"_"$sampleId"_QC.txt
echo -e "$rawSequenceQuality" >> "$seqId"_"$sampleId"_QC.txt

#Upload to BaseSpace
#TODO

#clean up
rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq.gz "$seqId"_"$sampleId"_"$laneId"_R2.fastq.gz