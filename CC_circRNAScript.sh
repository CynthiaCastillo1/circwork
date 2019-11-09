#!/bin/bash -l
#SBATCH --job-name=MycircRNAScript
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --workdir="/athena/elementolab/scratch/cyc4001/SLURM_Worklog"
#SBATCH --output=MycircRNAScript_%A_%a.out
#SBATCH --error=MycircRNAScript_%A_%a.err
#SBATCH --mem=50G

##SBATCH --mem-per-cpu=




spack load -r bzip2@1.0.6+shared%gcc@4.8.5
#spack load -r bzip2

output_dir=/athena/elementolab/scratch/cyc4001/DLBCL_analysis


path="/athena/elementolab/scratch/akv3001/DLBCL_Analysis/STRANDED_DLBCL_CellLines/FASTQ_DLBCLStranded_CellLines"
##path to all samples


echo "Path=$path"
echo $SLURM_TASK_ID

file=$(ls --ignore="*.gz" ${path}| tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Changed to SLURM  parallel processing environment variable

cd $TMPDIR  #copy the samples foldeR

echo "Processing $file"

Sample=$(basename "$file")

Sample=${Sample%.*}

echo "COPYING Sample Name = $Sample"

rsync -r -v $path/$file/*fastq.gz $TMPDIR

echo "Merging multi lane fastq if any.."

mkdir $TMPDIR/${Sample}_merged

      echo "Inside $file"

      for R in 1 2; do

              zcat $TMPDIR/*R$R*.fastq.gz | gzip -n -9 > $TMPDIR/${Sample}_merged/"$Sample"_R$R.fastq.gz &

              echo "..."    

	  done

        wait

rm $TMPDIR/*.gz

echo 'Processing' ${Sample}

        #--------Formatting multiple lane Data for Tophat---------

        F1="$TMPDIR/${Sample}_merged/"$Sample"*R1*.fastq.gz"

        F2="$TMPDIR/${Sample}_merged/"$Sample"*R2*.fastq.gz"

        echo "F1= " ${F1}

        echo "F2= " ${F2}

FASTA='/athena/elementolab/scratch/akv3001/GenomeReference_hg19/circRNA_hg19_annotations/hg19.fa'
REF_annotation='GenomeReference_hg19/circRNA_hg19_annotations/hg19_ref.txt'


##-------- Create storage directory on node--------------------------------------------
        mkdir $TMPDIR/${Sample}_STAR

#------------STAR Alignment Command--------------------------------#

STAR --chimSegmentMin 10 --runThreadN 10 \
     --genomeDir /athena/elementolab/scratch/akv3001/GenomeReference_hg19/Indexed_genome/hg19_Genomedir/ \
     --readFilesIn ${F1} ${F2} \
     --readFilesCommand zcat \
   ##--outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${Sample} \

ls -lh

CIRCexplorer2 parse $TMPDIR/${Sample}_STAR/ -t STAR *Chimeric.out.junction > CIRCexplorer2_parse.log

##CIRCexplorer2 parse -t STAR $TMPDIR/${Sample}_STAR/*Chimeric.out.junction > CIRCexplorer2_parse.log
##fast_circ.py parse -t STAR $TMPDIR/${Sample}_STAR/*Chimeric.out.junction > CIRCexplorer2_parse.log


##fast_circ.py annotate -r $REF_annotation -g $FASTA \
## -b ${Sample}_back_spliced_junction.bed \
## ${Sample}_CircRNA_Known.txt --low-confidence  > CIRCexplorer2_${Sample}.annotate.log
CIRCexplorer2 annotate -r $REF_annotation -g $FASTA \
  -b ${Sample}_back_spliced_junction.bed \
 ${Sample}_CircRNA_Known.txt --low-confidence  > CIRCexplorer2_${Sample}.annotate.log


mv low_conf_circularRNA_known.txt ${Sample}_low_conf_circularRNA_known.txt
##create file or get path?

#output_dir=$(dirname $PATH)  ##$output_dir
output_dir=/athena/elementolab/scratch/cyc4001
mkdir ${output_dir}/Circexplorer2_Summary
mkdir ${output_dir}/Circexplorer2_Summary/STAR_Bams
mkdir ${output_dir}/Circexplorer2_Summary/STAR_fusion
mkdir ${output_dir}/Circexplorer2_Summary/circRNA_Results
mkdir ${output_dir}/Circexplorer2_Summary/Result_Logs


##$TMPDIR/${Sample}=$PATH

echo "$TMPDIR" $TMPDIR

#rsync -v $TMPDIR/${Sample}/_STAR/*.bam ${output_dir}/Circexplorer2_Summary/STAR_Bams
#rsync -r -v $TMPDIR/${Sample}/_STAR/*Chimeric.out.junction  ${output_dir}/Circexplorer2_Summary/STAR_fusion
#rsync -r -v $TMPDIR/*.txt ${output_dir}/Circexplorer2_Summary/circRNA_Results
#rsync -r -v $TMPDIR/*.logs ${output_dir}/Circexplorer2_Summary/Result_Logs
#rsync -r -v $TMPDIR/denovo ${output_dir}/Circexplorer2_Summary/Denovo_Results/${Sample}
#rsync -r -v $TMPDIR/abs ${output_dir}/Circexplorer2_Summary/Denovo_Results/${Sample}
#rsync -r -v $TMPDIR/as ${output_dir}/Circexplorer2_Summary/Denovo_Results/${Sample}


