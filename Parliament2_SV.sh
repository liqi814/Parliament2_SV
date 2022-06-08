#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2022-04-20 11:44
# * Last modified : 2022-04-20 11:44
# * Filename      : Parliament2_SV.sh
# * Description   : fix docker container data volumn issue. best version

SAMPLE_ID=$1
mkdir $SAMPLE_ID

sudo docker container prune -f

cp /nfs/projects/dbGap/hg38_bam_file/${SAMPLE_ID}/${SAMPLE_ID}.merge.marked_duplicates.bam input_files/input.bam
cp /nfs/projects/dbGap/hg38_bam_file/${SAMPLE_ID}/${SAMPLE_ID}.merge.marked_duplicates.bai input_files/input.bam.bai

mkdir -p ${PWD}/mnt/docker
#sudo docker volume create  --opt type=none --opt o=bind --opt device=${PWD}/input_files in
sudo docker volume create  --opt type=none --opt o=bind --opt device=${PWD}/${SAMPLE_ID} out
sudo docker volume create  --opt type=none --opt o=bind --opt device=${PWD}/mnt/docker dnanexus

#sudo docker run -v ${PWD}/mnt/docker:/home/dnanexus -v ${PWD}/input_files/:/home/dnanexus/in/ -v ${PWD}/${SAMPLE_ID}/:/home/dnanexus/out/ dnanexus/parliament2:latest --bam input.bam --bai input.bam.bai  -r Homo_sapiens_assembly38.fasta  --fai Homo_sapiens_assembly38.fasta.fai --prefix  ${SAMPLE_ID}.${SubChrom} --filter_short_contigs --manta --genotype --delly_deletion --delly_insertion --delly_inversion --delly_duplication --lumpy

sudo docker run --mount source=dnanexus,target=/home/dnanexus/ --mount source=in,target=/home/dnanexus/in/ --mount source=out,target=/home/dnanexus/out/ dnanexus/parliament2:latest --bam input.bam --bai input.bam.bai -r Homo_sapiens_assembly38.fasta --fai Homo_sapiens_assembly38.fasta.fai --prefix ${SAMPLE_ID} --filter_short_contigs --manta --genotype --delly_deletion --delly_insertion --delly_inversion --delly_duplication --lumpy

#sudo docker run -v ${PWD}:/home/dnanexus   -it --entrypoint sh dnanexus/parliament2:latest
#sudo docker run -v ${PWD}:/home/dnanexus dnanexus/parliament2:latest rm -rf mnt/docker

# sudo docker container stop $CONTAINERID
sudo docker container prune -f
sudo docker volume rm dnanexus
sudo docker volume rm out
