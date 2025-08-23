


# sanyuan_gwas


dislocker-fuse -V /dev/sde2 --user-password="FSZsy09876" -- /mnt/bitlocker/
dislocker-fuse -V /dev/sdf2 --user-password="FSZsy09876" -- /mnt/bitlocker2/
mount -t ntfs-3g -o loop,rw /mnt/bitlocker/dislocker-file /home/guanguiwen/20T1
mount -t ntfs-3g -o loop,rw /mnt/bitlocker2/dislocker-file /home/guanguiwen/20T2
sanyuan
ls *_R1.clean.fastq.gz| awk -F'_' '{print $1" "$1"_"$2}' > sample

第一步 trim
thread=40
cat sample | while read line; do
    arr=($line)
    sample=${arr[0]}
    id=${arr[1]}
fastp \
-i ${id}_R1.fastq.gz \
-I ${id}_R2.fastq.gz \
-o /home/guanguiwen/data1/${id}_R1.clean.fastq.gz \
-O /home/guanguiwen/data1/${id}_R2.clean.fastq.gz \
--thread $thread \
--length_required 36 \
--html $id.fastp_report.html \
--json $id.fastp_report.json
done

第二步 bwa
ref=/home/guanguiwen/10T4/sanyuan_ref/REF/Homo_sapiens_assembly19.fasta
thread=46
output=/home/guanguiwen/10T5/sanyuan_bam/
while read -r sample id
do
RG='@RG\tID:'$id'\tSM:'$sample'\tLB:'$sample'\tPU:'$id'\tPL:illumina\tCN:genetron'
echo "Processing sample: $sample with ID: $id with RG:$RG"
bwa mem -R "$RG"   -t $thread        -M $ref      "${id}_R1.clean.fastq.gz"    "${id}_R2.clean.fastq.gz" | samtools view -bS - > ${output}${sample}.bam
done < sample
done


while read -r sample id; do RG='@RG\tID:'$id'\tSM:'$sample'\tLB:'$sample'\tPU:'$id'\tPL:illumina\tCN:genetron'; echo "Processing sample: $sample with ID: $id with RG:$RG"; bwa mem -R "$RG"   -t $thread        -M $ref      "${id}_R1.clean.fastq.gz"    "${id}_R2.clean.fastq.gz" | samtools view -bS - > ${output}${sample}.bam; done < sample

bwa mem -R "@RG\tID:K155378N-K155378N_E250066377-L01\tSM:K155378N\tLB:K155378N\tPU:K155378N_E250066377-L01\tPL:illumina\tCN:genetron" \
          -t 8 -M /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R1.fq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R2.fq.gz \
          >/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.bam 

第三步 排序
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /picard_tools/SortSam.jar \
          I=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.bam \
          O=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          MAX_RECORDS_IN_RAM=5000000 SO=coordinate VALIDATION_STRINGENCY=SILENT \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G4_bam_sort/G4_bam_sort-k155378n-node4.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam

第四步 interval
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T RealignerTargetCreator \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list.log \
          -known /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -known /idc_obu/pipeline/Database/reference_hg19/Mills_and_1000G_gold_standard.indels.b37.vcf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G5_bam_interval/G5_bam_interval-k155378n-node5.sh.log 2>&1

第五步 markdup
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /picard_tools/MarkDuplicates.jar \
          INPUT=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          OUTPUT=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup.bam \
          METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G6_bam_markdup/G6_bam_markdup-k155378n-node6.sh.log 2>&1

第六步 indelrealign
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T IndelRealigner \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup.bam \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -targetIntervals /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam.log \
          -known /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -known /idc_obu/pipeline/Database/reference_hg19/Mills_and_1000G_gold_standard.indels.b37.vcf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G7_bam_indelrealign/G7_bam_indelrealign-k155378n-node7.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam

第七步 baserecal
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T BaseRecalibrator \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -knownSites /idc_obu/pipeline/Database/reference_hg19/dbsnp_138.b37.vcf \
          -knownSites /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp.log \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G8_bam_baserecal/G8_bam_baserecal-k155378n-node8.sh.log 2>&1

第八步 covariates
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -T AnalyzeCovariates \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -BQSR /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -plots /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp.pdf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G9_bam_covariates/G9_bam_covariates-k155378n-node9.sh.log 2>&1

第九步 printreads-recal
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -T PrintReads \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -BQSR /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned_recal.bam \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G10_bam_printreads/G10_bam_printreads-k155378n-node10.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned_recal.bam






第一步 trim
java -Xmx16g -jar /software/trimmomatic-0.33.jar PE -phred33 -trimlog \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/trim.log \
          -threads 8 \
		  /RawData_classed_H3C/Taq_MGI_T7/KY_WGS_120/K155378N_E250066377-L01_R1.fastq.gz \
          /RawData_classed_H3C/Taq_MGI_T7/KY_WGS_120/K155378N_E250066377-L01_R2.fastq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R1.fq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_unpaired_R1.fq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R2.fq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_unpaired_R2.fq.gz \
          ILLUMINACLIP:/data/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G2_trim/G2_trim-k155378n-node2.sh.log 2>&1

第二步 bwa
bwa mem -R "@RG\tID:K155378N-K155378N_E250066377-L01\tSM:K155378N\tLB:K155378N\tPU:K155378N_E250066377-L01\tPL:illumina\tCN:genetron" \
          -t 8 -M /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R1.fq.gz \
          /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/00_QC/K155378N_E250066377-L01_clean_R2.fq.gz \
          >/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.bam 

第三步 排序
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /picard_tools/SortSam.jar \
          I=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.bam \
          O=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          MAX_RECORDS_IN_RAM=5000000 SO=coordinate VALIDATION_STRINGENCY=SILENT \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G4_bam_sort/G4_bam_sort-k155378n-node4.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam

第四步 interval
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T RealignerTargetCreator \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list.log \
          -known /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -known /idc_obu/pipeline/Database/reference_hg19/Mills_and_1000G_gold_standard.indels.b37.vcf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G5_bam_interval/G5_bam_interval-k155378n-node5.sh.log 2>&1

第五步 markdup
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /picard_tools/MarkDuplicates.jar \
          INPUT=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_E250066377-L01.sort.bam \
          OUTPUT=/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup.bam \
          METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G6_bam_markdup/G6_bam_markdup-k155378n-node6.sh.log 2>&1

第六步 indelrealign
java -Xmx16g -Djava.io.tmpdir=/idc_obu/pipeline/WGS/normal_only/tmp -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T IndelRealigner \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup.bam \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -targetIntervals /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_intervals.list \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam.log \
          -known /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -known /idc_obu/pipeline/Database/reference_hg19/Mills_and_1000G_gold_standard.indels.b37.vcf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G7_bam_indelrealign/G7_bam_indelrealign-k155378n-node7.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam

第七步 baserecal
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -T BaseRecalibrator \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -knownSites /idc_obu/pipeline/Database/reference_hg19/dbsnp_138.b37.vcf \
          -knownSites /idc_obu/pipeline/Database/reference_hg19/1000G_phase1.indels.b37.vcf \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -log /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp.log \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G8_bam_baserecal/G8_bam_baserecal-k155378n-node8.sh.log 2>&1

第八步 covariates
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -T AnalyzeCovariates \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -BQSR /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -plots /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp.pdf \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G9_bam_covariates/G9_bam_covariates-k155378n-node9.sh.log 2>&1

第九步 printreads-recal
java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
          -T PrintReads \
          -R /idc_obu/pipeline/Database/reference_hg19/REF/Homo_sapiens_assembly19.fasta \
          -I /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned.bam \
          -BQSR /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_recal.grp \
          -o /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned_recal.bam \
          >>/idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/shell/G10_bam_printreads/G10_bam_printreads-k155378n-node10.sh.log 2>&1

samtools index /idc_obu/project/CHB/03.analysis/batch4/20250510/obu/K155378N/01_aln/K155378N_rmdup_realigned_recal.bam
