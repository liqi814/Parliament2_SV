#! /bin/bash
SAMPLE_ID=$1
REF_FILE="/nfs/projects/dbGap/hg38_bam_file/reference/Homo_sapiens_assembly38.fasta"

#https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
rm -rf $SAMPLE_ID;
mkdir $SAMPLE_ID;
cellflow_jobID=""
for FRIST_READ in  `ls /nfs/archive/p2018/FASTQ/GENOME/${SAMPLE_ID}/*/${SAMPLE_ID}*_R1_001.fastq.gz`
do
  LANE=`echo $FRIST_READ | grep -o 'L00[0-9]' | sed 's/L00//'`
  FLOWCELL_BARCODE=`echo $FRIST_READ | rev | cut -d '/' -f2 | rev`
  SECOND_READ=`echo $FRIST_READ | sed 's/_R1_/_R2_/'`
  RGINFO=`samtools view -H /nfs/archive/p2018/ALIGNMENT/BUILD37/DRAGEN/GENOME_AS_FAKE_EXOME/$SAMPLE_ID*/*.realn.recal.bam | grep '^@RG' | grep ${FLOWCELL_BARCODE}.${LANE}`
  LB=`echo $RGINFO | cut -d ' ' -f3`
  jobID=`echo "/nfs/goldstein/software/bwa/bwa mem $REF_FILE \
  -t 7 $FRIST_READ $SECOND_READ  \
  -R '@RG\tID:${FLOWCELL_BARCODE}.${LANE}\t${LB}\tPL:ILLUMINA\tPU:${FLOWCELL_BARCODE}.${LANE}\tSM:$SAMPLE_ID' | \
  samtools sort -@8 -o ${SAMPLE_ID}/${SAMPLE_ID}.${FLOWCELL_BARCODE}.${LANE}.bam -" |  \
  qsub -terse -cwd -V -S /bin/bash -j y -N ${SAMPLE_ID}.${FLOWCELL_BARCODE}.${LANE} -pe threaded  10`
  cellflow_jobID="${cellflow_jobID},${jobID}"
done

cellflow_jobID=`echo $cellflow_jobID | sed 's/^,//' | sed 's/,$//g'`

# merge BAM files
merge_jobID=`echo "samtools merge -@8 \
${SAMPLE_ID}/${SAMPLE_ID}.merge.bam \
${SAMPLE_ID}/${SAMPLE_ID}.*bam ; \
if [ -f \"${SAMPLE_ID}/${SAMPLE_ID}.merge.bam\" ]; then \
  ls ${SAMPLE_ID}/${SAMPLE_ID}.*bam | grep -v merge | xargs rm ; \
fi; " | \
qsub -terse -hold_jid ${cellflow_jobID} -cwd -V -S /bin/bash -j y -N ${SAMPLE_ID}_mergeBam  -pe threaded  8`

# Mark duplicates
echo "/home/am5153/bin/rtg-tools-3.9.1/jre/bin/java -jar /nfs/goldstein/software/picard-tools-2.23.8/picard.jar MarkDuplicates \
I=${SAMPLE_ID}/${SAMPLE_ID}.merge.bam \
O=${SAMPLE_ID}/${SAMPLE_ID}.merge.marked_duplicates.bam \
M=${SAMPLE_ID}/${SAMPLE_ID}.marked_dup_metrics.txt \
CREATE_INDEX=true ; \
if [ -f \"${SAMPLE_ID}/${SAMPLE_ID}.merge.marked_duplicates.bam\" ]; then \
  ls ${SAMPLE_ID}/${SAMPLE_ID}.merge.bam |  xargs rm ; \
fi; " | \
qsub -hold_jid ${merge_jobID} -cwd -V -S /bin/bash -j y -N ${SAMPLE_ID}_MarkDuplicates  -pe threaded 2
