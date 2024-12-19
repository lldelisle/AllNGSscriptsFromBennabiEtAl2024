# For cellplex:

# toolshed.g2.bx.psu.edu/repos/iuc/cite_seq_count/cite_seq_count/1.4.4+galaxy0
# command_version: CITE-seq-Count v1.4.4
CITE-seq-Count --threads ${GALAXY_SLOTS:-4} --read1 'input_R1.fastqsanger.gz' --read2 'input_R2.fastqsanger.gz' --tags 'CMO_seq csv.csv' \
    --cell_barcode_first_base 1 --cell_barcode_last_base 16 --umi_first_base 17 --umi_last_base 28 --bc_collapsing_dist 1 \
    --umi_collapsing_dist 2  --expected_cells 24000 --whitelist 'cellranger_barcodes_3M-february-2018.txt.txt' --max-error 2    
 gunzip Results/read_count/barcodes.tsv.gz 
 gunzip Results/read_count/features.tsv.gz 
 gunzip Results/read_count/matrix.mtx.gz 
 gunzip Results/umi_count/barcodes.tsv.gz 
 gunzip Results/umi_count/features.tsv.gz 
 gunzip Results/umi_count/matrix.mtx.gz

# For Gene expression:

# toolshed.g2.bx.psu.edu/repos/iuc/rna_starsolo/rna_starsolo/2.7.10b+galaxy3
STAR  --runThreadN ${GALAXY_SLOTS:-4} --genomeLoad NoSharedMemory --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' \
    --sjdbOverhang 100 --sjdbGTFfile 'input.gtf'   --soloType CB_UMI_Simple   --readFilesIn input_R1.fastq.gz input_R2.fastq.gz --soloCBmatchWLtype 1MM_multi  \
    --readFilesCommand zcat   --soloCBwhitelist 'cellranger_barcodes_3M-february-2018.txt.txt' --soloBarcodeReadLength 1 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 \
    --soloUMIlen 12   --soloStrand Forward --soloFeatures Gene --soloUMIdedup 1MM_CR --soloUMIfiltering - --quantMode TranscriptomeSAM GeneCounts  \
    --outSAMattributes NH HI AS nM GX GN CB UB --outSAMtype BAM SortedByCoordinate  --soloCellFilter None  --soloOutFormatFeaturesGeneField3 'Gene Expression' \
    --outSAMunmapped None --outSAMmapqUnique 60  --limitOutSJoneRead 1000 --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000     
 mv Solo.out/Gene Solo.out/soloFeatures 
 cat <(echo "Barcodes:") Solo.out/Barcodes.stats <(echo "Genes:") Solo.out/soloFeatures/Features.stats > 'RNA STARSolo on data 47, data 38, and others: Barcode/Feature Statistic Summaries.txt'  
 samtools view -b -o 'RNA STARSolo on data 47, data 38, and others: Alignments.bam' Aligned.sortedByCoord.out.bam

# toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/dropletutils/1.10.0+galaxy2
# command_version:
mkdir 'tenx.input' 
 ln -s 'RNA STARSolo on data 47, data 38, and others: Matrix Gene Counts raw.mtx' 'tenx.input/matrix.mtx' 
 ln -s 'RNA STARSolo on data 47, data 38, and others: Genes raw.tsv' 'tenx.input/genes.tsv' 
 ln -s 'RNA STARSolo on data 47, data 38, and others: Barcodes raw.tsv' 'tenx.input/barcodes.tsv' 
  mkdir 'tenx.output' 
  Rscript '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/a9caad671439/dropletutils/scripts/dropletutils.Rscript' '/data/galaxy/galaxy/jobs/000/107/107308/configs/tmptrjig6gs'

