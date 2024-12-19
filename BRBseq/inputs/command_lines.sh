# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2
# command_version:GNU Awk 4.2.0, API: 2.0
# awk program is  $1 != "Name"{print $2} 
env -i $(which awk) --sandbox -v FS='	' -v OFS='	' --re-interval -f "/data/galaxy/galaxy/jobs/000/102/102043/configs/tmp5ftor_rg" "barcodes_96_V5A_brb.txt.tabular" > "Text reformatting on data 4.tabular"

# toolshed.g2.bx.psu.edu/repos/iuc/rna_starsolo/rna_starsolo/2.7.10b+galaxy1
# command_version:
STAR  --runThreadN ${GALAXY_SLOTS:-4} --genomeLoad NoSharedMemory --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 100 --sjdbGTFfile 'mm10_custom102_allGastruloids_min10_extended.gtf.gtf'   --soloType CB_UMI_Simple   --readFilesIn Lib_A_S1_R2.fastq.gz.fastqsanger.gz Lib_A_S1_R1.fastq.gz.fastqsanger.gz --soloCBmatchWLtype 1MM_multi  --readFilesCommand zcat   --soloCBwhitelist 'Text reformatting on data 4.tabular' --soloBarcodeReadLength 1 --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 14 --soloAdapterSequence '-' --soloAdapterMismatchesNmax 1 --clipAdapterType Hamming   --soloStrand Forward --soloFeatures Gene --soloUMIdedup 1MM_All --quantMode TranscriptomeSAM GeneCounts --outSAMattributes NH HI AS nM CR UR GX GN CB UB sM sS sQ --outSAMtype BAM SortedByCoordinate  --soloCellFilter None  --soloOutFormatFeaturesGeneField3 'Gene Expression'  '--outSAMunmapped Within' --outSAMmapqUnique 60  --limitOutSJoneRead 1000 --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000     
 mv Solo.out/Gene Solo.out/soloFeatures 
 cat <(echo "Barcodes:") Solo.out/Barcodes.stats <(echo "Genes:") Solo.out/soloFeatures/Features.stats > 'RNA STARSolo on data 7, data 2, and others: Barcode/Feature Statistic Summaries.txt'   
 samtools view -b -o 'RNA STARSolo on data 7, data 2, and others: Alignments.bam' Aligned.sortedByCoord.out.bam
