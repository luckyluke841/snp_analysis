#!/usr/bin/env bash

#  processzips.sh

#  Reads must be included as _R1 and _R2
#  See loopfiles.sh and email_loopfiles for multiple samples.

#################################################################################
#  Dependencies ---
#    bwa mem
#    samtools
#    abyss
#    gatk
#    pilon
#    igvtools
#    gff if available for annotation
#    bamtools
#    python 3, modules: sys, csv, xlsxwriter
#    Bruc_MLST.sh
#    spoligoSpacerFinder.sh
#    trimming sequence: /usr/local/bin/bbmap/resources/nextera.fa.gz
#   File containing high quality SNPs, Volumes/Mycobacterium/Go_To_File/HighestQualitySNPs.vcf
#   Reference in fasta format, /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta
#################################################################################

###################
# HELP
###################
function help () {

printf "\n\n*** Help ***\n"
printf "See inside script for list of dependencies\n\n"
printf "Usage: \n"
printf "\tworking directory containing zipped paired FASTQ files\n"
printf "\tDo not provide a sample_type if more than 1 sample is in working directory\n"
printf "\tOptions:\n"
printf "\t\t-h help\n"
printf "\t\t-d debug\n"
printf "\t\t-e email\n"
printf "\t\t\tcalling email will loop multiple samples\n"
printf "\t\t\tcall with all, tod, jess, suelee, or off\n"
printf "\t\t-t sample_type\n"
printf "\t\t\tsample_type options:\n"
printf "\t\t\t\tTBBOV, H37Rv\n"
printf "\t\t\t\tab1, mel1, mel2, mel3, suisall, suis1, suis2, suis3, suis4, suis5, canis, ceti1, ceti2, ovis\n"
printf "\t\t\t\tpara, past, h5n2 secd, taylorella\n\n"
printf "example:\n"
printf "\t~$ processZips.sh -t TBBOV \t\t# working directory with one sample\n"
printf "\t~$ processZips.sh -e tod \t\t# loop and send email to tod.p.stuber@usda.gov\n"
printf "\t~$ processZips.sh -e all \t\t# loop and send email to all\n\n"

exit 1
}

###################
# RUN_SAMPLE
###################
function run_sample() {

alias pause='read -p "$LINENO Enter"'
echo "current directory"
pwd
startingdir=`pwd`

# Move zip files to their own directory
mkdir ./zips
mv *.fastq* ./zips
mkdir ./alignment
cd alignment/

# Make alias links in BWA-GATK directory to zip files
ls ../zips/*.fastq* | while read file; do ln -s $file; done
pause

if [ $sample_type == ab1 ]; then
    cp /home/shared/brucella/abortus1/script_dependents/NC_00693c.fasta ./
    hqs="/home/shared/brucella/abortus1/script_dependents/NC_00693cHighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/abortus1/newFiles"
    sharedSAN="/home/shared/brucella/abortus1/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == ab3 ]; then
    cp /home/shared/brucella/abortus3/script_dependents/CP007682-7683c.fasta ./
    hqs="/home/shared/brucella/abortus3/script_dependents/CP007682-7683cHighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/abortus3/newFiles"
    sharedSAN="/home/shared/brucella/abortus3/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################
    
elif [ $sample_type == mel1 ]; then
    cp /home/shared/brucella/melitensis1/script_dependents/mel-bv1-NC003317.fasta ./
    hqs="/home/shared/brucella/melitensis1/script_dependents/mel-bv1-NC003317-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/melitensis-bv1/newFiles"
#sharedSAN="/home/shared/brucella/melitensis/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == mel2 ]; then
    cp /home/shared/brucella/melitensis2/script_dependents/mel-bv2-NC012441.fasta ./
    hqs="/home/shared/brucella/melitensis2/script_dependents/mel-bv2-NC012441-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/melitensis-bv2/newFiles"
#sharedSAN="/home/shared/brucella/melitensis/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == mel3 ]; then
    cp /home/shared/brucella/melitensis3/script_dependents/mel-bv3-NZCP007760.fasta ./
    hqs="/home/shared/brucella/melitensis3/script_dependents/mel-bv3-NZCP007760-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/melitensis-bv3/newFiles"
#sharedSAN="/home/shared/brucella/melitensis/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suisall ]; then
    cp /home/shared/brucella/suis5/script_dependents/NZ_CP00771c.fasta ./
    hqs="/home/shared/brucella/suis5/script_dependents/B-513-highqualitysnps.vcf"
    #bioinfo="/bioinfo11/TStuber/Results/brucella/suisall/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suis1 ]; then
    cp /home/shared/brucella/suis1/script_dependents/NC_01725c.fasta ./
    hqs="/home/shared/brucella/suis1/script_dependents/NC_01725cHighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/suis1/newFiles"
    sharedSAN="/home/shared/brucella/suis1/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suis2 ]; then
    cp /home/shared/brucella/suis2/script_dependents/Bsuisbv2-94-11.fasta ./
    hqs="/home/shared/brucella/suis2/script_dependents/suis2HighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/suis2/newFiles"
    sharedSAN="/home/shared/brucella/suis2/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suis3 ]; then
    cp /home/shared/brucella/suis3/script_dependents/B-REF-BS3-686.fasta ./
    hqs="/home/shared/brucella/suis3/script_dependents/B15-0007-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/suis3/newFiles"
    sharedSAN="/home/shared/brucella/suis3/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suis4 ]; then
    cp /home/shared/brucella/suis4/script_dependents/B-REF-BS4-40.fasta ./
    hqs="/home/shared/brucella/suis4/script_dependents/suis4HighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/suis4/newFiles"
    sharedSAN="/home/shared/brucella/suis4/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == suis5 ]; then
    cp /home/shared/brucella/suis5/script_dependents/NZ_CP00771c.fasta ./
    hqs="/home/shared/brucella/suis5/script_dependents/B-513-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/suis5/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == canis ]; then
    cp /home/shared/brucella/canis/script_dependents/BcanisATCC23365.fasta ./
    hqs="/home/shared/brucella/canis/script_dependents/canisHighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/canis/newFiles"
    sharedSAN="/home/shared/brucella/canis/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == ceti1 ]; then
    cp /home/shared/brucella/ceti1/script_dependents/Bceti1Cudo.fasta ./
    hqs="/home/shared/brucella/ceti1/script_dependents/ceti1HighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/ceti1/newFiles"
    sharedSAN="/home/shared/brucella/ceti1/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == ceti2 ]; then
    cp /home/shared/brucella/ceti2/script_dependents/Bceti2-TE10759.fasta ./
    hqs="/home/shared/brucella/ceti2/script_dependents/ceti2HighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/ceti2/newFiles"
    sharedSAN="/home/shared/brucella/ceti2/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == ovis ]; then
    cp /home/shared/brucella/ovis/script_dependents/BovisATCC25840.fasta ./
    hqs="/home/shared/brucella/ovis/script_dependents/BovisATCC25840HighestQualitySNPs.vcf"
    bioinfo="/bioinfo11/TStuber/Results/brucella/ovis/newFiles"
    sharedSAN="/home/shared/brucella/ovis/newFiles"

    # Run BrucMLST.sh
    echo "Starting Bruc_MLST.sh"
    cd ../zips
    ${BRUC_MLST} &
    cd ../alignment/
    echo "Moving forward from Bruc_MLST.sh"

    ###################################################################

elif [ $sample_type == taylorella ]; then
    cp /home/shared/taylorella/NC018108.fasta ./
    hqs="/home/shared/taylorella/TE-004-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/gen-bact/taylorella/newFiles"

elif [ $sample_type == tay1 ]; then
    cp /home/shared/taylorella/09-0932.fasta ./
    hqs="/home/shared/taylorella/TE-004-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/gen-bact/taylorella/newFiles"

elif [ $sample_type == tay5 ]; then
    cp /home/shared/taylorella/92-0972-DFS5.fasta ./
    hqs="/home/shared/taylorella/15-0094-taylorella-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/gen-bact/taylorella/newFiles"

    ###################################################################

    # Lineage Bov-Afri, AF2122
elif [ $sample_type == TBBOV ]; then
    cp /home/shared/mycobacterium/tbc/snppipeline/tbbov/NC_002945.fasta ./
    hqs="/home/shared/mycobacterium/tbc/snppipeline/tbbov/HighestQualitySNPs.vcf"
    gff_file="/home/shared/mycobacterium/tbc/snppipeline/tbbov/NC_002945.gff"
    bioinfo="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles"
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

# Lineage Bov-Afri, AF2122-updated_2017-01-25
elif [ $sample_type == AF2122 ]; then
    cp /home/shared/mycobacterium/tbc/snppipeline/tbbov/NC_002945.fasta ./
    hqs="/home/shared/mycobacterium/tbc/snppipeline/tbbov/HighestQualitySNPs.vcf"
    gff_file="/home/shared/mycobacterium/tbc/snppipeline/tbbov/NC_002945.gff"
    bioinfo="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles"
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################

    # H37Rv general
elif [ $sample_type == H37Rv ]; then
    cp /home/shared/mycobacterium/tbc/snppipeline/tb4b/NC000962.fasta ./
    hqs="/home/shared/mycobacterium/tbc/snppipeline/tb4b/15-3162-highqualitysnps.vcf"
    bioinfo="/bioinfo11/TStuber/Results/mycobacterium/tbc/h37/newfiles"
    gff_file="/home/shared/mycobacterium/tbc/snppipeline/tb4b/NC_000962.gff"

    # Run spoligoSpacerFinder.sh
    echo "Starting spoligoSpacerFinder.sh"
    ${SPOLIGOSPACERFINDER} &
    echo "Moving forward from spoligoSpacerFinder.sh"

    ###################################################################
    ###################################################################

elif [ $sample_type == para ]; then
    cp /home/shared/mycobacterium/mott/paratb/NC_002944.fasta ./
    hqs="/home/shared/mycobacterium/mott/paratb/HQ-NC002944.vcf"
    bioinfo="/bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/newFiles"
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

elif [ $sample_type == past ]; then
    cp /home/shared/pasteurella/NZ_CM001581.fasta ./
    hqs="/home/shared/pasteurella/BTYP-9814-pasteurella-highqualitysnps.vcf"
    #bioinfo="/bioinfo11/TStuber/Results/gen-bact/Pasteurella/newFiles"
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

elif [ $sample_type == h5n2 ]; then
    cp /home/shared/virus/ai/h5n2/TY-BC-FAV10-2014.fasta ./
    hqs="/home/shared/virus/ai/h5n2/14111-1-highqualitysnps.vcf"
    bioinfo="/bioinfo11/MKillian/Analysis/results/influenza/h5n2/snp_analysis/newfiles/"
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

elif [ $sample_type == h5nx ]; then
    cp /home/shared/virus/ai/h5nx/BC-turkey-PB2-HA-MP.fasta ./
    hqs="/home/shared/virus/ai/h5nx/11602-1-highqualitysnps.vcf"
    bioinfo=""
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

    ###################################################################

elif [ $sample_type == secd ]; then
    cp /home/shared/virus/secd/scriptDependents/KC210145.fasta ./
    hqs="/home/shared/virus/secd/scriptDependents/HighestQualitySNPs.vcf"
    #bioinfo=""
    #sharedSAN="/home/shared/mycobacterium/bovis/newFiles"

else
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel1, mel2, mel3, suisall, suis1, suis2, suis3, suis4, suis5, canis, ceti1, ceti2, ovis, TBBOV, H37Rv, para, past, h5n2 secd, taylorella"
    exit 1
fi
pause

echo "*** Trimming"

forReads=`ls | grep _R1`
echo "Forward Reads to be trimmed: $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads to be trimmed:: $revReads"

#Trim the reads with bbmap tool kit (bbduk plugin)
#about twice as fast as trimmomatic

strain=$(echo $revReads | sed 's/_.*//' | sed 's/\..*//')
echo -e "Quality trimming sample "$strain""

    ${BBDUK} -Xmx80g \
    in1="$forReads" \
    in2="$revReads" \
    ref="/usr/local/bin/bbmap/resources/nextera.fa.gz" \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=lr trimq=5 \
    minlen=36 \
    out1=trimmed_reads/${strain}_Trimmed_R1.fastq.gz \
    out2=trimmed_reads/${strain}_Trimmed_R2.fastq.gz \
    stats=trim_stats.txt \
    qchist=qc_by_base.txt \
    threads=auto \
    showspeed=f

mv -v trimmed_reads/${strain}_Trimmed_R1.fastq.gz ./
mv -v trimmed_reads/${strain}_Trimmed_R2.fastq.gz ./
rm -r trimmed_reads
rm "$forReads"
rm "$revReads"

forReads=`ls | grep _R1`
echo "Forward Reads to be used after trimmed: $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads to be used after trimmed:: $revReads"

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

#   retrieves reference name and name from sorted BAM file name
r=`echo $ref | sed 's/\..*//'`
n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

echo "***Reference naming convention:  $r"
echo "***Isolate naming convention:  $n"

${SAMTOOLS} faidx $ref
java -Xmx4g -jar ${PICARD} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
    echo "Index and dict are present, continue script"
    else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    ${SAMTOOLS} faidx $ref
    java -Xmx4g -jar ${PICARD} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
        if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
        fi
fi

# See echo comments
echo "***bwa index $r"
${BWA} index $ref

# -t sets the number of threads/cores
# -r ST	 Specify the read group in a format like ‘@RG\tID:foo\tSM:bar’ Needed for GATK
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
echo "***Making Sam file"
${BWA} mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam

# -b	 Output in the BAM format.
# -h	 Include the header in the output.
#-F INT	 Skip alignments with bits present in INT [0]
echo "***Making Bam file"
${SAMTOOLS} view -bh -F4 -T $ref $n.sam > $n.raw.bam

if [ $sample_type == secd ]; then
        echo "secd, not assembling unmapped reads"
else
####### unmapped reads #######
#Bam with mapped and unmapped reads
${SAMTOOLS} view -bh -T $ref $n.sam > $n.all.bam
#Strip off the unmapped reads
${SAMTOOLS} view -h -f4 $n.all.bam > $n.unmappedReads.sam
#Create fastqs of unmapped reads to assemble
java -Xmx4g -jar ${PICARD} SamToFastq INPUT=$n.unmappedReads.sam FASTQ=${n}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-unmapped_R2.fastq
rm $n.all.bam
rm $n.unmappedReads.sam
${ABYSS} name=${n}_abyss k=64 in="${n}-unmapped_R1.fastq ${n}-unmapped_R2.fastq"

mkdir ../unmappedReads
mv ${n}-unmapped_R1.fastq ../unmappedReads
mv ${n}-unmapped_R2.fastq ../unmappedReads
mv ${n}_abyss-3.fa ../unmappedReads
mv ${n}_abyss-8.fa ../unmappedReads
mv ${n}_abyss-stats ../unmappedReads
mv *coverage* ../unmappedReads
rm *abyss*
######################
fi

echo "***Sorting Bam"
${SAMTOOLS} sort $n.raw.bam -o $n.sorted.bam
echo "***Indexing Bam"
${SAMTOOLS} index $n.sorted.bam
# Remove duplicate molecules

echo "***Marking Duplicates"
java -Xmx4g -jar  ${PICARD} MarkDuplicates INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $n.dup.bam"
${SAMTOOLS} index $n.dup.bam

# Run Pilon
mkdir pilon
java -Xmx16G -jar ${PILON} --genome $ref --bam ${n}.dup.bam --output ./pilon/${n}-pilon --vcf --vcfqe --tracks --iupac
awk ' $5 != "." || $7 != "PASS" {print $0}' ./pilon/${n}-pilon.vcf | awk '$7 != "LowCov" {print $0}' | awk '$5 != "<DUP>" {print $0}' | awk '$8 !~ /SVTYPE=INS/ {print $0}' | awk '$8 !~ /SVTYPE=DEL/ {print $0}' | awk '$8 !~ /SVTYPE=DUP/ {print $0}' > ${n}.pilon.vcf
cp ${n}.pilon.vcf ./pilon

##########

#if [ $gff_file ]; then
#    echo "Annotating $n.pilon.vcf"
#    awk '$3 == "gene" {print $0}' $gff_file | awk '{print $4, $5, $9}' > list.genes
#
#    while read l; do
#    echo $l | awk '{for(i=$1;i<=$2;i++) print i, $3}'
#    done < list.genes | sed -e 's/\([0-9]*\).*;\(Name=.*\);gbkey.*\(gene_biotype=.*\);\(locus_tag=.*\)/\1   \2;\3;\4/' > expand.gene
#
#    #Split header lines from position calls
#    grep "#" $n.pilon.vcf > header
#    grep -v "#" $n.pilon.vcf > body
#
#    #http://unix.stackexchange.com/questions/113898/how-to-merge-two-files-based-on-the-matching-of-two-columns
#
#    awk 'BEGIN{OFS="\t"}NR==FNR {h[$1] = $2; next} {print $1,$2,h[$2],$4,$5,$6,$7,$8,$9,$10}' expand.gene body | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' > body2
#
#    cat header body2 > ${n}.pilon-annotated.vcf
#
#    rm body
#    rm body2
#    rm header
#    rm expand.gene
#    rm list.genes
#
#else
#    echo "gff file not given"
#fi

echo "***Deleting Files"
rm $n.sam
rm $n.raw.bam
rm $n.dup.bam
rm $n.dup.bam.bai
rm $n.sorted.bam
rm $n.sorted.bam.bai
rm $forReads
rm $revReads
rm igv.log

rm $n.SNPsMapzeroNoHeader.vcf
rm $n.vcfFormated
rm $n.keepTheseZeroCovPositions
rm $n.duplicates
rm $n.zeroCoveragePositions
rm $n.calledSNPpositions
rm $n.body
rm $n.header
rm $n.hapreadyOnlySNPs.vcf

###################################
# The next 5 steps collect metrics
###################################

##Quality Score Distribution
#echo "***Quality Score Distribution"
#java -Xmx4g -jar ${PICARD} QualityScoreDistribution REFERENCE_SEQUENCE=$ref INPUT=$n.dup.bam CHART_OUTPUT=$n.QualityScorceDistribution.pdf OUTPUT=$n.QualityScoreDistribution ASSUME_SORTED=true
#
##Mean Quality by Cycle
#echo "***Mean Quality by Cycle"
#java -Xmx4g -jar ${PICARD} CollectMultipleMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.dup.bam OUTPUT=$n.Quality_by_cycle PROGRAM=MeanQualityByCycle ASSUME_SORTED=true
#
##Collect Alignment Summary Metrics
#echo "***Collect Alignment Summary Metrics"
#java -Xmx4g -jar ${PICARD} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.dup.bam OUTPUT=$n.AlignmentMetrics ASSUME_SORTED=true
#
##Collect GC Bias Error
#echo "***Collect GC Bias Error"
#java -Xmx4g -jar ${PICARD} CollectGcBiasMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.dup.bam OUTPUT=$n.CollectGcBiasMetrics CHART_OUTPUT=$n.GC.PDF ASSUME_SORTED=true
#
##Collect Insert Size Metrics
#echo "***Collect Insert Size Metrics"
#java -Xmx4g -jar ${PICARD} CollectInsertSizeMetrics REFERENCE_SEQUENCE=$ref INPUT=$n.dup.bam HISTOGRAM_FILE=$n.InsertSize.pdf OUTPUT=$n.CollectInsertSizeMetrics ASSUME_SORTED=true
#
#cat $n.AlignmentMetrics >> $n.Metrics_summary.xls
#cat $n.CollectInsertSizeMetrics >> $n.Metrics_summary.xls
#
#echo "***Organizing files"
#
##Move to qualityvalues subfolder
#mkdir qualityvalues
#mv $n.GC.PDF qualityvalues/
#mv $n.QualityScorceDistribution.pdf qualityvalues/
#mv $n.InsertSize.pdf qualityvalues/
#mv $n.Quality_by_cycle*.pdf qualityvalues/

#Remove files
rm $n.CollectInsertSizeMetrics
rm $n.Quality_by_cycle.quality_distribution_metrics
rm $n.Quality_by_cycle.quality_by_cycle_metrics
rm $n.Quality_by_cycle.alignment_summary_metrics
rm $n.CollectGcBiasMetrics
rm $n.QualityScoreDistribution
rm $n.coverage

###########################
echo "***Getting stats for $n"

echo "fastq.gz file sizes:" > $n.stats2.txt
ls -lh ../zips/*gz | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt
read1_size=`ls -lh ../zips/*gz | awk '{print $5}' | egrep -v '^$' | head -1`
read2_size=`ls -lh ../zips/*gz | awk '{print $5}' | egrep -v '^$' | tail -1`

echo "Unmapped fastq.gz file sizes:" >> $n.stats2.txt
gzip ../unmappedReads/*fastq
ls -lh ../unmappedReads/*gz | awk '{print $5}' | egrep -v '^$' >> $n.stats2.txt

echo "Unmapped contig count:" >> $n.stats2.txt
grep -c ">" ../unmappedReads/${n}_abyss-3.fa >> $n.stats2.txt
unmapped_contig_count=`grep -c ">" ../unmappedReads/${n}_abyss-3.fa`

echo "" >> $n.stats2.txt
sed -n 7,8p $n.FilteredReads.xls | awk '{print $2}' >> $n.stats2.txt
unpaired_reads=`sed -n 7,8p $n.FilteredReads.xls | awk '{print $2}' | tail -1`
sed -n 7,8p $n.FilteredReads.xls | awk '{print $3}' >> $n.stats2.txt
total_reads_pairs=`sed -n 7,8p $n.FilteredReads.xls | awk '{print $3}' | tail -1`
sed -n 7,8p $n.FilteredReads.xls | awk '{print $8}' >> $n.stats2.txt
duplicate_reads=`sed -n 7,8p $n.FilteredReads.xls | awk '{print $8}' | tail -1`
readcount=`sed -n 8p $n.FilteredReads.xls | awk '{print $3}'`
echo "" >> $n.stats2.txt

echo "***Bamtools running"
aveCoverage=`${BAMTOOLS} coverage -in $n.dup.bam | awk '{sum+=$3} END { print sum/NR"X"}'`
aveCoveragenoX=`echo $aveCoverage | sed 's/X//'`
echo "Average depth of coverage: $aveCoverage" >> $n.stats2.txt

#genome coverage
percGenomeMissing=`awk -v x="$zeroposition" -v y="$refsize" 'BEGIN { print(x/y)*100}'`
percGenomeCoverage="$(echo "100 - $percGenomeMissing" | bc)"
echo "Percent of reference with coverage: ${percGenomeCoverage}%" >> $n.stats2.txt

#cat $n.stats2.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" >> $n.stats.txt
cat $n.stats2.txt > $n.stats.txt
echo "" >> $n.stats.txt

rm $n.stats2.txt
###########################

#  Add Insert_Size and Read_Length to stats.txt file
echo 'Mean_Insert_Size  Standard_Deviation:' >> $n.stats.txt
awk 'BEGIN {OFS="\t"} { print $5,$6 }' $n.Quality_by_cycle.insert_size_metrics | awk 'FNR == 8 {print $0}' >> $n.stats.txt

echo 'Mean_Read_Length:' >> $n.stats.txt
awk 'BEGIN {OFS="\t"} { print $16 }' $n.AlignmentMetrics | awk 'FNR == 10 {print $0}' >> $n.stats.txt
average_read_length=`awk 'BEGIN {OFS="\t"} { print $16 }' $n.AlignmentMetrics | awk 'FNR == 10 {print $0}' | sed 's/\..*//'`
echo "" >> $n.stats.txt

#  Add SNP call numbers to stats.txt file
echo "SNP and zero coverage positions in $n.pilon.vcf:" >> $n.stats.txt
egrep -v "#" $n.pilon.vcf | grep -c ".*" >> $n.stats.txt

echo "SNPs of AC2 and QUAL > 300:" >> $n.stats.txt
egrep -v "#" $n.pilon.vcf | egrep "AC=2" | awk '$6 > 300' | grep -c ".*" >> $n.stats.txt
quality_snps=`egrep -v "#" $n.pilon.vcf | egrep "AC=2" | awk '$6 > 300' | grep -c ".*"`

#  Show Mean Coverage at Terminal and coverageReport
echo "Mean Coverage"

echo "Sample identified and ran as:  $sample_type" >> $email_summary_top

echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% "
echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% " >> /scratch/report/coverageReport.txt
echo -e "$n \t $readcount \t ${aveCoverage} \t ${percGenomeCoverage}% " >> $email_summary_top

mv $n.Metrics_summary.xls qualityvalues/
mv $n.stats.txt qualityvalues/
mv $n.FilteredReads.xls qualityvalues/
rm $n.Quality_by_cycle.insert_size_metrics
rm $n.AlignmentMetrics
mv ${startingdir}/fastq ${startingdir}/spoligo
rm ${startingdir}/spoligo/*fastq
rm -r ${startingdir}/temp
ln qualityvalues/$n.stats.txt ./stats-$n.txt

cp $0 ./

#echo "***Sending files to the Network"
#cp -r ${startingdir} ${bioinfo}

printf "%s\t%s\t%s\t%s\t%'d\t%'d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\n" $sample_type $n $read1_size $read2_size $total_reads_pairs $unpaired_reads $duplicate_reads $average_read_length $r $aveCoveragenoX $percGenomeCoverage $unmapped_contig_count $quality_snps >> /scratch/report/pre_stat_table.txt
printf "%s\t%s\t%s\t%s\t%'d\t%'d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\n" $sample_type $n $read1_size $read2_size $total_reads_pairs $unpaired_reads $duplicate_reads $average_read_length $r $aveCoverage $percGenomeCoverage $unmapped_contig_count $quality_snps >> /scratch/report/stat_table_cumulative.txt

#Make dailyStats.txt for each stats.txt made for each isolate.
echo "" >> $email_summary_bottom
echo "" >> $email_summary_bottom
echo "" >> $email_summary_bottom
echo "ADD_MARKER" >> $email_summary_bottom
echo "" >> $email_summary_bottom
echo "<------- $n $sample_type ------->" >> $email_summary_bottom
cat qualityvalues/$n.stats.txt >> $email_summary_bottom
cat qualityvalues/$n.stats.txt >> /home/shared/stats

echo "**************************** END $n ****************************"

#
#  Created by Stuber, Tod P - APHIS on 11/08/12.
#

}

###################
# OLIGO_IDENTIFIER
###################
function oligo_identifier () {

sample_name=`echo $1 | sed 's:[./]::g'`
oligo_path=`pwd`
log_oligo="${oligo_path}/log_oligo_${sample_name}.txt"

onemismatch () {
patt=($1)
for ((i=0; i<${#patt[0]}; i++)); do
    patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
regex=$(IFS='|'; echo "${patt[*]}")

}

#################################################################################

echo "**********************ID FILES START**********************"

# updated oligos, 2014-12-11
pab1="AATTGTCGGATAGCCTGGCGATAACGACGC"
pab3="CACACGCGGGCCGGAACTGCCGCAAATGAC"
pab5="GCTGAAGCGGCAGACCGGCAGAACGAATAT"
pmel="TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
ether2="CGAAATCGTGGTGAAGGACGGGACCGAACC" # added 2017-01-13, 0 means mel bv3
p63B1="CCTGTTTAAAAGAATCGTCGGAACCGCTCT" # added 2017-01-13, 0 means mel bv2
p16M0="TCCCGCCGCCATGCCGCCGAAAGTCGCCGT" # added 2017-01-13, 0 means mel bv1
psuis1="TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
psuis2="GGCAATCATGCGCAGGGCTTTGCATTCGTC"
psuis3="CAAGGCAGATGCACATAATCCGGCGACCCG"
pceti1="GTGAATATAGGGTGAATTGATCTTCAGCCG"
pceti2="TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
pcanis4="CTGCTACATAAAGCACCCGGCGACCGAGTT"
pcanis="ATCGTTTTGCGGCATATCGCTGACCACAGC"
povis="CACTCAATCTTCTCTACGGGCGTGGTATCC"

# updated oligos, 2015-02-02

primertb157="CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
primertb7="TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
primertbbov="CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
primertb5="CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
primertb2="ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
primertb3="GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
primertb4="CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
primertb6="ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"

# paraTB specific
ppara="CCTTTCTTGAAGGGTGTTCG|CGAACACCCTTCAAGAAAGG"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

n=`echo $forReads | sed 's/_.*//' | sed 's/\..*//'`
echo "working on: $n"

ab1=`grep -c $pab1 $forReads`
ab3=`grep -c $pab3 $forReads`
ab5=`grep -c $pab5 $forReads`
mel=`grep -c $pmel $forReads`
count_ether2=`grep -c $ether2 $forReads`
count_p63B1=`grep -c $p63B1 $forReads`
count_p16M0=`grep -c $p16M0 $forReads`
suis1=`grep -c $psuis1 $forReads`
suis2=`grep -c $psuis2 $forReads`
suis3=`grep -c $psuis3 $forReads`
ceti1=`grep -c $pceti1 $forReads`
ceti2=`grep -c $pceti2 $forReads`
canis4=`grep -c $pcanis4 $forReads`
canis=`grep -c $pcanis $forReads`
ovis=`grep -c $povis $forReads`
para=`egrep -ic $ppara $forReads`

onemismatch $primertb157
tb157=`egrep -c $regex $forReads`
echo "tb157 $tb157"

onemismatch $primertb7
tb7=`egrep -c $regex $forReads`
echo "tb7 $tb7"

onemismatch $primertbbov
tbbov=`egrep -c $regex $forReads`
echo "tbbov $tbbov"

onemismatch $primertb5
tb5=`egrep -c $regex $forReads`
echo "tb5 $tb5"

onemismatch $primertb2
tb2=`egrep -c $regex $forReads`
echo "tb2 $tb2"

onemismatch $primertb3
tb3=`egrep -c $regex $forReads`
echo "tb3 $tb3"

onemismatch $primertb4
tb4=`egrep -c $regex $forReads`
echo "tb4 $tb4"

onemismatch $primertb6
tb6=`egrep -c $regex $forReads`
echo "tb6 $tb6"

bruccounts=`echo "$ab1 $ab3 $ab5 $mel $suis1 $suis2 $suis3 $ceti1 $ceti2 $canis4 $canis $ovis $count_ether2 $count_p63B1 $count_p16M0"`
tbcounts=`echo "$tb157 $tb7 $tbbov $tb5 $tb2 $tb3 $tb4 $tb6"`
paracounts=`echo "$para"`
echo $bruccounts
echo $tbcounts
echo $paracounts

brucbinary=`echo $bruccounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
tbbinary=`echo $tbcounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
parabinary=`echo $paracounts | awk '{for(i=1;i<=NF;i++) if ($i >= 10) print 1; else print 0}' | tr -cd "[:print:]"`

echo $brucbinary
echo $tbbinary
echo $parabinary

#count every occurance of 1 in binary.
check=`echo $brucbinary | grep -o "1" | wc -l`
echo "Brucella check= $check"

if [[ $check > 0 ]]; then
	echo "Brucella species found"
	tagname=`grep $n /bioinfo11/TStuber/Results/brucella/bruc_tags.txt`
	i=$brucbinary

	if [ $i == 111111111111111 ]
        then
		echo "*** Odd Isolate, Unexpected findings ***" >> $log_oligo
	    exit 1

	elif [ $i == 011111111111111 ]
        then
		echo "Brucella abortus bv 1, 2 or 4" >> $log_oligo
        sample_type="ab1"

	elif [ $i == 101111111111111 ]
		then
		echo "Brucella abortus bv 3" >> $log_oligo
        sample_type="ab3"

	elif [ $i == 110111111111111 ]
    	then
        echo "Brucella abortus bv 5, 6 or 9" >> $log_oligo
        sample_type="ab1"

	elif [ $i == 111011111111110 ] || [ $i == 000001010110110 ]
        then
        echo "Brucella melitensis" >> $log_oligo
        sample_type="mel1"

	elif [ $i == 111011111111101 ] || [ $i == 000001010110100 ] || [ $i == 010001010110100 ] || [ $i == 111001111110101 ]
        then
        echo "Brucella melitensis" >> $log_oligo
        sample_type="mel2"

	elif [ $i == 111011111111011 ] || [ $i == 111001111110011 ]
        then
        echo "Brucella melitensis" >> $log_oligo
        sample_type="mel3"

	elif [ $i == 111101111111111 ]
        then
		echo "Brucella suis bv1" >> $log_oligo
        sample_type="suis1"

	elif [ $i == 111110111111111 ]
        then
		echo "Brucella suis bv2" >> $log_oligo
        sample_type="suis2"

	elif [ $i == 111111011111110 ]
        then
		echo "Brucella suis bv3" >> $log_oligo
        sample_type="suis3"

	elif [ $i == 111111101111111 ] || [ $i == 111111100111111 ]
    	then
        echo "Brucella ceti 1" >> $log_oligo
        sample_type="ceti1"

	elif [ $i == 111111110111111 ]
    	then
        echo "Brucella ceti 2" >> $log_oligo
        sample_type="ceti2"

	elif [ $i == 111111111011110 ]
    	then
        echo "Brucella suis bv4" >> $log_oligo
        sample_type="suis4"

	elif [ $i == 111111111001110 ]
    	then
        echo "Brucella canis" >> $log_oligo
        sample_type="canis"

	elif [ $i == 111111111110111 ]
    	then
        echo "Brucella ovis" >> $log_oligo
        sample_type="ovis"

	else
		echo "*** Odd Isolate, Unexpected findings, See /home/shared/brucella/bruc_oligo_identifier_output.txt ***"
    	echo "***bruc_oligo_identifier cannot find a pattern for $n, see line $LINENO of script***"
		echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts}, Binary: ${brucbinary} ${tbbinary}" >> $email_summary_top
	fi
fi

#count every occurance of 1 in binary.
check=`echo $tbbinary | grep -o "1" | wc -l`
echo "TB complex check= $check"

if [[ $check > 0 ]]; then
        echo "TB complex species found"
        tagname=`grep $n /bioinfo11/TStuber/Results/mycobacterium/Untitled.txt`
	i=$tbbinary

        if [ $i == 11101111 ] || [ $i == 11101101 ]
        then
        sample_type="H37Rv"
        echo "TB1" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo
	
        elif [ $i == 01100111 ]
        then
        sample_type="H37Rv"
        echo "TB2" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo 

        elif [ $i == 01101011 ] || [ $i == 11101011 ]
        then
        sample_type="H37Rv"
        echo "TB3" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo

        elif [ $i == 01101111 ]
        then
        sample_type="H37Rv"
        echo "TB4a" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo


        elif [ $i == 01101101 ] || [ $i == 11101101 ] || [ $i == 01101111 ]
        then
        sample_type="H37Rv"
        echo "TB4b" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        "$tbcounts" >> $log_oligo

        elif [ $i == 11111111 ]
        then
        sample_type="H37Rv"
        echo "TB5" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo

        elif [ $i == 11001111 ]
        then
        sample_type="H37Rv"
        echo "TB6" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo

        elif [ $i == 10101110 ]
        then
        sample_type="H37Rv"
        echo "TB7 are as TBBOV" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo

        elif [ $i == 11001110 ] || [ $i == 11011110 ] || [ $i == 11001100 ]
        then
        sample_type="TBBOV"
        echo "TBBOV" >> $log_oligo
        echo "$tbbinary" >> $log_oligo
        echo "$tbcounts" >> $log_oligo

        else
        sample_type="No match found"
        echo "***oligo_identifier cannot place $n with TB reference, see line $LINENO of script***"
		echo "Oligo counts:  ${tbcounts}, Binary:  $tbbinary"
		echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts}, Binary: ${brucbinary} ${tbbinary}" >> $email_summary_top
        exit 1
	fi
fi

#count every occurance of 1 in binary.
check=`echo $parabinary | grep -o "1" | wc -l`
echo "M. paratb check= $check"

if [[ $check > 0 ]]; then
echo "MAC species found"
tagname=`grep $n /bioinfo11/TStuber/Results/mycobacterium/Untitled.txt`
i=$parabinary

    if [ $i == 1 ]; then
        sample_type="para"
        echo "M. paratuberculosis found"
        echo "M. paratuberculosis" >> $log_oligo
        echo "$parabinary" >> $log_oligo
        echo "$paracounts" >> $log_oligo
    else
        echo "oligo_identifier.sh could not find a match for $n"
        echo "oligo_identifier.sh could not find a match for $n" >> $email_summary_top
        echo "${n} Unable to find a reference, oligo_identifier.sh stats: Oligo counts: ${bruccounts} ${tbcounts} ${paracounts}, Binary: ${brucbinary} ${tbbinary} ${parabinary}" >> $email_summary_top
    fi
fi

allbinary=`echo ${brucbinary} ${tbbinary} ${parabinary} | grep -c "1"`
if [[ $allbinary == 0  ]]; then 
	echo "${n} NEEDS SPECIAL ATTENTION!!!" >> $email_summary_top
	echo "PLEASE GIVE TOD SPECIAL INSTRUCTIONS FOR ${n}.  This sample was NOT identified as Brucella, TB complex or avium complex.  If you know something about this isolate please send me an email.  CC all on email.  Thanks!" | mutt -s "${n} NEEDS SPECIAL ATTENTION!!!" -- "tod.p.stuber@usda.gov suelee.robbe-austerman@aphis.usda.gov" #Doris.M.Bravo@aphis.usda.gov john.b.fevold@aphis.usda.gov Patrick.M.Camp@aphis.usda.gov David.T.Farrell@aphis.usda.gov" 
fi

echo "Sample ${n}, ${tagname}, Oligo counts: Bruc ${bruccounts} TB ${tbcounts} MAC ${paracounts}, Binary: Bruc ${brucbinary} TB ${tbbinary} MAC ${parabinary}, ID:  ${sample_type}"

#Push to logfile
echo "Sample ${n}, ${tagname}, Oligo counts: Bruc ${bruccounts} TB ${tbcounts} MAC ${paracounts}, Binary: Bruc ${brucbinary} TB ${tbbinary} MAC ${parabinary}, ID:  ${sample_type}" >> $log_oligo

mv $log_oligo  ../
cd ..

#cannot excute alias inside function
#pause
#read -p "$LINENO Enter"

run_sample $sample_type 2>&1 | tee log_processZips_${sample_name}.txt

}

# END OF FUCTIONS

#################################################################################
######                                                                     ######
######                            # START #                                ######
######                                                                     ######
#################################################################################

echo "**************************************************************************"
echo "**************************** START ${PWD##*/} ****************************"
echo "**************************************************************************"

PICARD=`which picard.jar`
if [[ -z $PICARD ]]; then
    echo "picard.jar not in PATH"
    echo "picard version >1.14"
    echo "Add picard.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

BWA=`which bwa`
if [[ -z $BWA ]]; then
    echo "bwa is not in PATH"
    echo "Add bwa to PATH"
    echo "See line: $LINENO"
    exit 1
fi

SAMTOOLS=`which samtools`
    if [[ -z $SAMTOOLS ]]; then
    echo "samtools is not in PATH"
    echo "Add samtools to PATH"
    echo "See line: $LINENO"
    exit 1
fi

ABYSS=`which abyss-pe`
    if [[ -z $ABYSS ]]; then
    echo "abyss-pe is not in PATH"
    echo "Add abyss-pe to PATH"
    echo "See line: $LINENO"
    exit 1
fi

BAMTOOLS=`which bamtools`
if [[ -z $BAMTOOLS ]]; then
    echo "Bamtools is not in PATH"
    echo "Add Bamtools to PATH"
    echo "See line: $LINENO"
    exit 1
fi

GATK=`which GenomeAnalysisTK.jar`
if [[ -z $GATK ]]; then
    echo "GenomeAnalysisTK.jar is not in PATH"
    echo "Add GenomeAnalysisTK.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

IGVTOOLS=`which igvtools.jar`
if [[ -z $IGVTOOLS ]]; then
    echo "igvtools.jar is not in PATH"
    echo "Add igvtools.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

PILON=`which pilon.jar`
if [[ -z $PILON ]]; then
    echo "pilon.jar is not in PATH"
    echo "Add pilon.jar to PATH"
    echo "See line: $LINENO"
exit 1
fi

BBDUK=`which bbduk.sh`
if [[ -z $BBDUK ]]; then
    echo "igvtools.jar is not in PATH"
    echo "Add igvtools.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

SPOLIGOSPACERFINDER=`which spoligoSpacerFinder.sh`
if [[ -z $SPOLIGOSPACERFINDER ]]; then
    echo "spoligoSpacerFinder.sh is not in PATH"
    echo "Add spoligoSpacerFinder.sh to PATH"
    echo "See line: $LINENO"
    exit 1
fi

BRUC_MLST=`which Bruc_MLST.sh`
if [[ -z $BRUC_MLST ]]; then
    echo "spoligoSpacerFinder.sh is not in PATH"
    echo "Add spoligoSpacerFinder.sh to PATH"
    echo "See line: $LINENO"
    exit 1
fi


email_summary_top="/scratch/report/dailyReport.txt"
# Detail daily stats shown in email (bottom portion of email)
email_summary_bottom="/scratch/report/dailyStats.txt"

root=`pwd`

shopt -s expand_aliases
alias pause='read -p "$LINENO Enter"'

# Set flags

email=
sample_type=
debug=
while getopts ':hde:t:' OPTION; do
    case $OPTION in
        h) hflag=1
        ;;
        d) debug=1
        ;;
        e) email=$OPTARG
        ;;
        t) sample_type=$OPTARG
        ;;
        ?) printf "\n\nIncorrect argument given"; help
        ;; 
    esac
done
shift $(($OPTIND - 1))

if [[ $hflag ]]; then
    help
    exit 1
fi

printf "\nARE YOU RUNNING AS ROOT?\n\n"
echo "hflag: $hflag"
echo "debug: $debug"
echo "email: $email"
printf "sample_type: $sample_type\n\n"
#read -p "*** Check Arguments and Press Enter to Continue ***"

if [[ $sample_type ]]; then
    count=`ls | grep -c "_R1"`
    if [[ $count -gt 1 ]]; then
        debug=1
        printf "\n\nEXPECTING ONLY A SINGLE SAMPLE\n"
        printf "ONLY A SINGLE SAMPLE IS ALLOWED WHEN USING -t OPTION\n\n"
        printf "\nExited at $LINENO\n"
        help
        exit 1
    fi
fi

##############################
for i in *.fastq.gz; do
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    echo "n is : $n"
    mkdir -p $n
    mv $i $n/
done

if [[ $email ]] && [[ $sample_type ]]; then
    printf "\n\nSAMPLE_TYPE CANNOT BE SPECIFIED WHEN LOOPING SAMPLES BY CALLING EMAIL\n"
    printf "SAMPLE_TYPE MUST BE DETERMINED BY OLIGO_IDENTIFIER FUCTION\n"
    printf "\nExited at $LINENO\n"
    help
    exit 1
fi

if [[ $sample_type ]]; then
    cd ${root}/$n
    run_sample $sample_type | tee tee_processZips_out.txt
fi

# If email then loop samples
# Or if single and no sample_type given
if [[ $email ]] || [[ -z $sample_type ]]; then

currentdir=`pwd`

printf "ref_type\tsample\tR1_zip\tR2_zip\ttotal_read_prs\tup_reads\t%%dup_reads\tave_read_length\tref\tave_cov_X\tper_cov\tunmapped_contigs\tquality_snps\n" > /scratch/report/stat_table.txt
printf "" > /scratch/report/pre_stat_table.txt
printf "\n" >> /scratch/report/stat_table_cumulative.txt
date >> /scratch/report/stat_table_cumulative.txt

echo "Start Time: `date`" > /scratch/report/dailyTime
starttime=`date +%s`

echo "Please wait.  Searching for TB complex, Brucella and paratuberculosis oligos and then starting appropriate processZips.sh argument"

#`loopFiles.sh` &&
date >> /scratch/report/coverageReport.txt

# Reset spoligo and bruc mlst check file
echo "" > /scratch/report/spoligoCheck.txt
echo "" > /scratch/report/mlstCheck.txt
echo "WG Spoligo Check" >> /scratch/report/spoligoCheck.txt
echo "Brucella MLST Check" >> /scratch/report/mlstCheck.txt

#Reset file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "TB Number" "Octal Code" > "/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

echo "Isolate Total_Bases AveDep %>Q15" | awk '{ printf("%-12s %-12s %-10s %-10s\n", $1, $2, $3, $4) }' > $email_summary_top
printf "\n\n" >> $email_summary_top

echo ""  > $email_summary_bottom

for sample_directory in ./*/; do
	echo "The cd is $currentdir"
	cd $currentdir
	echo "$sample_directory started"
	cd ./$sample_directory
	mkdir ./temp
	cp *R1*.fastq.gz ./temp
    cd ./temp
	gunzip *.fastq.gz

    if [[ $debug ]]; then
    # do not put into a subprocess, iterate one at a time
        oligo_identifier $sample_directory
    else
        oligo_identifier $sample_directory &
    fi
done
wait
###

echo "" >> /scratch/report/dailyTime
echo "End Time: `date`" >> /scratch/report/dailyTime
endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> /scratch/report/dailyTime

echo "e-mailing files"

cat /scratch/report/dailyTime > /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt
echo "ADD_MARKER" >> /scratch/report/email_processZips.txt
cat $email_summary_top >> /scratch/report/email_processZips.txt

cat /scratch/report/spoligoCheck.txt >> /scratch/report/email_processZips.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/email_processZips.txt

echo "ADD_MARKER" >> /scratch/report/email_processZips.txt

cat $email_summary_bottom >> /scratch/report/email_processZips.txt
echo "" >> /scratch/report/email_processZips.txt

grep -v '*' /scratch/report/email_processZips.txt | grep -v "Stats for BAM file" | sed 's/ADD_MARKER/******************************************/g' > /scratch/report/email_processZips2.txt

# Create "here-document"
cat >${root}/excelwriterstats.py <<'EOL'
#!/usr/bin/env python

import sys
import csv
import xlsxwriter

filename = sys.argv[1].replace(".txt",".xlsx")
wb = xlsxwriter.Workbook(filename)
ws = wb.add_worksheet("Sheet1")
with open(sys.argv[1],'r') as csvfile:
    table = csv.reader(csvfile, delimiter='\t')
    i = 0
    for row in table:
        ws.write_row(i, 0, row)
        i += 1

col = len(row)
print (filename, ":", i, "x", col)

wb.close()

EOL

chmod 755 ${root}/excelwriterstats.py

sort -k1,2 /scratch/report/pre_stat_table.txt >> /scratch/report/stat_table.txt

${root}/excelwriterstats.py /scratch/report/stat_table.txt

column -t /scratch/report/stat_table.txt > /scratch/report/stat_table.temp; mv /scratch/report/stat_table.temp /scratch/report/stat_table.txt
enscript /scratch/report/stat_table.txt -B -j -r -f "Courier7" -o - | ps2pdf - /scratch/report/stat_table.pdf

if [[ $email -eq 1 ]]; then
	email_list="tod.p.stuber@usda.gov"
elif [[ $email == "all" ]]; then
    email_list="tod.p.stuber@usda.gov"
elif [[ $email == "tod" ]]; then
    email_list="tod.p.stuber@usda.gov"
elif [[ $email == "jess" ]]; then
    email_list="tod.p.stuber@usda.gov"
elif [[ $email == "suelee" ]]; then
    email_list="tod.p.stuber@usda.gov"
elif [[ $email == "off" ]]; then
    email_list="off"
    printf "\n\n\tNo email sent\n\n"
else
    echo "Email was not sent: email argument: $email"
fi

cat /scratch/report/email_processZips2.txt | mutt -a /scratch/report/stat_table.xlsx /scratch/report/stat_table.pdf -s "WGS results" -- $email_list

date >> /scratch/report/mlstCheck_all.txt
cat /scratch/report/mlstCheck.txt >> /scratch/report/mlstCheck_all.txt

rm ${root}/excelwriterstats.py


fi

##############################










