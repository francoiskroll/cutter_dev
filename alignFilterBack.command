#!/bin/bash

### arguments are:
    # -c path to config.csv, c for config
    # -r path to folder with fastq reads, r for reads
    # -a path to folder with fasta reference sequences, a for alignment
    # -u if present, we only have one fastq file per sample, rather than one with Fw reads and one with Rv reads
    # -l whether to filter the bam file or not; l for fiLter; filtering is ON if flag is given, OFF if not

# in addition, arguments from filterBAM.command:
    # -i = input = bam file to process
    # -e = PhrEd score = minimum Phred score
    # -f = floor = minimum read span
    # -s = maximum proportion soft-clipped
    # -d = padding around double-strand break site
    # -p = keep only primary alignments? yes or no
    # -o = output = name of bam file to output

shopt -s nullglob

######################################

### function checkPath to check paths given
checkPath() {
    local file_path="$1"

    # Check if the file exists
    if [ -e "$file_path" ]; then
        :
    else
        echo "Error. File does not exist: $file_path"
        exit 1
    fi
}

### function checkDir to check path leads to an existing directory
checkDir() {
    local folder_path="$1"

    # Check if the file exists
    if [ -d "$folder_path" ]; then
        :
    else
        echo "Error. Folder does not exist: $folder_path"
        exit 1
    fi
}

######################################

### read the flags/arguments given by user
# each : means it requires an argument
# -u and -l do not require one, so no : after
while getopts c:r:a:ule:f:s:d:p: flag
do
    case "${flag}" in
        c) configpath=${OPTARG};;
        r) fastqpath=${OPTARG};;
        a) refpath=${OPTARG};;
        u) unique=${OPTARG};;
        l) filter=${OPTARG};;
        e) min_phred=${OPTARG};;
        f) min_readspan=${OPTARG};;
        s) max_softprop=${OPTARG};;
        d) dsbpad=${OPTARG};;
        p) primary=${OPTARG};;
    esac
done

### check paths given by user
checkPath "$configpath"
checkDir "$fastqpath"
checkDir "$refpath"

# if all OK, give a summary to user:
echo ""
echo "> Config file: $configpath"
echo "> Reads fastq: $fastqpath"
echo "> References fasta: $refpath"

# first convert config.xlsx to temporary config.csv
# note: reading a CSV file I created myself seems impossible
# if want to switch to using a CSV directly, still use ssconvert command below (from CSV to CSV)
ssconvert "$configpath" config.csv

################################################################################
### about directories

### we set main output folder as folder where the config file is
outdir=$(dirname "$configpath")

### create new folder for bam files we create
# will put folder in folder which contains the config file
mkdir "$(echo "$outdir"$"/bam")"

### if filtering is ON, create new folder for filtered BAMs & for new back-converted fastq (g for gunzip)
if [[ $* == *-l* ]]
then
    mkdir "$(echo "$outdir"$"/bamfilt")"
    mkdir "$(echo "$outdir"$"/filterfastqg")"
fi

################################################################################
### loop through rows of the config file        
while IFS=',' read -r well ref
do
    ### skip header row
    if [ "$well" == "well" ]; then
        continue
    fi

    ### tell user what we read
    echo ""
    echo "Well = $well"
    echo "Ref = $ref"

    ### in fastq folder, find files that start with well
    reads=$(find "$fastqpath" -type f -name "$well*")

    # check that exactly two fastq files were found
    # except if 'unique' flag is on (means we have only one fastq file per sample)
    readscount=$(echo "$reads" | wc -l)

    ### if unique file... ###
    if [[ $* == *-u* ]]
    then
        if [ "$readscount" -ne 1 ]; then
        echo "Error: -u flag is present so there should be exactly ONE fastq files starting with '$well', but found $readscount."
        exit 1
        fi
        
        # we have fastq file already (as there is only one), variable $reads
        # keep only filename
        firstnm=$(basename "$reads")

        echo "Forward reads = $reads"

        ### find fasta reference in the folder given by user
        refp=$(find "$refpath" -type f -name "$ref") # reference path

        ### in summary:
        echo "Aligning $reads to $refp"
        echo ""

        ### index the fasta reference
        bwa index -a bwtsw "$refp"

        ### do the alignment using bwa
        # create the output bam filename, which should be well_ref.bam
        # ref from config file should be ampliconname.fa
        # from there, will get ampliconname to append to output filenames
        amp="$(echo "$ref" | cut -d'.' -f 1)" # cut gets everything before first '.', so e.g. psen1_2.fa becomes psen1_2
        # paste well and amp together, so we have e.g. A01_psen1_2
        bam="$(echo "$well"$"_""$amp"$".bam")"
        # we will put the bam file in the new folder we created above called 'bam'
        bamp="$(echo "$outdir"$"/bam/""$bam")"

        # will also sort with samtools
        bwa mem -t 16 "$refp" "$reads" \
        | samtools sort > "$bamp"

    ### if fw/rv files... ###
    else # i.e. -u is not present, so we expect two fastq per sample
        if [ "$readscount" -ne 2 ]; then
        echo "Error: there should be exactly TWO fastq files starting with '$well', but found $readscount."
        exit 1
        fi

        # will assume one is _R1 (or _1), one is _R2 (or _2)
        # but will not assume which is which
        first=$(echo "$reads" | head -n 1)
        second=$(echo "$reads" | head -n 2 | tail -n 1)
        # keep only filename
        firstnm=$(basename "$first")
        # if first file has _R1 or _1 in its name
        if [[ "$firstnm" == *"_R1"* || "$firstnm" == *"_1"* ]]; then
            fwd="$first"
            rvs="$second"
        else
            rvs="$first"
            fwd="$second"
        fi

        echo "Forward reads = $fwd"
        echo "Reverse reads = $rvs"

        ### find fasta reference in the folder given by user
        refp=$(find "$refpath" -type f -name "$ref") # reference path

        ### in summary:
        echo "Aligning $fwd & $rvs to $refp"
        echo ""

        ### index the fasta reference
        bwa index -a bwtsw "$refp"

        ### do the alignment using bwa
        # create the output bam filename, which should be well_ref.bam
        # ref from config file should be ampliconname.fa
        # from there, will get ampliconname to append to output filenames
        amp="$(echo "$ref" | cut -d'.' -f 1)" # cut gets everything before first '.', so e.g. psen1_2.fa becomes psen1_2
        # paste well and amp together, so we have e.g. A01_psen1_2
        bam="$(echo "$well"$"_""$amp"$".bam")"
        # we will put the bam file in the new folder we created above called 'bam'
        bamp="$(echo "$outdir"$"/bam/""$bam")"

        # will also sort with samtools
        bwa mem -t 16 "$refp" "$fwd" "$rvs" \
        | samtools sort > "$bamp"
    fi

    # check we did create a bam file
    checkPath "$bamp"

    # index the sorted bam
    samtools index "$bamp"

    ################################
    ### filter the bam file

    # we filter if flag -l was given
    if [[ $* == *-l* ]]
    then
    
        # filtered bam to write is:
        # well_referenceamplicon_filt.bam
        bamfilt="$(echo "$well"$"_""$amp"$"_filt.bam")"
        # we will put the filtered bam file in a new folder we created above called 'bamfilt'
        bamfiltp="$(echo "$outdir"$"/bamfilt/""$bamfilt")"

        echo ""
        echo ">>> Filtering "$bam" into "$bamfilt""

        # filtering uses filterBam.command
        filterBam.command -i "$bamp" \
            -a "$refpath" \
            -e "$min_phred" \
            -f "$min_readspan" \
            -s "$max_softprop" \
            -d "$dsbpad" \
            -p "$primary" \
            -o "$bamfiltp"

        # filterBam.command includes sorting and indexing
        # check that we did create the filtered BAM file
        checkPath "$bamfiltp"

        ################################
        ### convert back to fastq
        # only do this if filtering is ON

        # we want Forward reads to be well_referenceamplicon_R1.fastq.gz, e.g. A01_psen1_2_R1.fa
        fwdfq="$(echo "$well"$"_""$amp"$"_R1.fastq")" # name of the file without .gz
        # put it in new folder filterfastqg
        fwdfqp="$(echo "$outdir"$"/filterfastqg/""$fwdfq")"

        # we want Reverse reads to be well_referenceamplicon_R2.fastq.gz, e.g. A01_psen1_2_R2.fa
        rvsfq="$(echo "$well"$"_""$amp"$"_R2.fastq")" # name of the file without .gz
        rvsfqp="$(echo "$outdir"$"/filterfastqg/""$rvsfq")"

        echo
        echo
        echo "---- [ CONVERTING BACK TO FASTQ"$" $bamfilt ] ----"
        echo "---- [ INTO"$" $fwdfq ] ----"
        echo
        echo

        # filterBam already sorts by chromosome coordinates
        # but for bamtofastq below it is important to sort by readnames
        # this is so each read is followed by its mate (i.e. forward read is followed by its reverse)
        # prepare the file name
        bamfiltNS="$(echo "$well"$"_""$amp"$"_tmp.bam")" # this is BAM from above, name-sorted
        # full path is
        bamfiltNSp="$(echo "$outdir"$"/bamfilt/""$bamfiltNS")"
        # now sort using samtools
        samtools sort -n "$bamfiltp" -o "$bamfiltNSp"

        # ready to convert
        # ! if unique -u, below will throw every read because it cannot find mates
        # use slightly different command
        if [[ $* == *-u* ]]
        then
            samtools fastq -0 "$fwdfqp" "$bamfiltNSp" # bedtools bamtofastq seems unable to deal with this case!
            # it doubles the read count even when outputting to a single fastq file
            checkPath "$fwdfqp" # check that we did create filtered Forward file
            # now gunzip them
            # it will delete original
            gzip -f "$fwdfqp"

        else
            bedtools bamtofastq -i "$bamfiltNSp" \
                -fq  "$fwdfqp" -fq2 "$rvsfqp"
            checkPath "$fwdfqp" # check that we did create filtered Forward file
            checkPath "$rvsfqp" # check that we did create filtered Reverse file
            # now gunzip them
            # it will delete original
            gzip -f "$fwdfqp"
            gzip -f "$rvsfqp"
        fi

        # remove temporary file
        rm "$bamfiltNSp"

    ### if no filtering, then we are done
    else
        echo ""
        echo ">>> Filtering is OFF."

        exit

    fi

    ################################

done < "config.csv"

### remove config.csv
rm config.csv

### copy folder filterfastqg to filterfastq
cp -rf "$(echo "$outdir"$"/filterfastqg")" "$(echo "$outdir"$"/filterfastq")"
# then gunzip all the files in filterfastq
gunzip "$(echo "$outdir"$"/filterfastq")"/*.gz

shopt -u nullglob
exit