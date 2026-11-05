> conda env list

##m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq.gz在/staging/biology/u6807805/Tala的資料夾
> zcat m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq.gz | head
> conda activate seqkit
> seqkit stats m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq.gz
##result##
file                                                      format  type  num_seqs        sum_len  min_len   avg_len  max_len
m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq.gz  FASTQ   DNA    209,331  3,335,501,175      382  15,934.1   44,552

> /staging/biology/u6807805/Tala2
> conda activate smudgeplot
> gunzip /staging/biology/u6807805/Tala/m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq.gz
>FastK -T32 -v -t1 -M120 -k21 m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq -NFastk_Table_bc1002

##################Hifiasm genome assembly################
> vim hifiasm.sh
-------------_______貼上以下內容________--------------
#!/bin/bash -l
#SBATCH -A MST109178
#SBATCH -p ngs372G
#SBATCH -c 56
#SBATCH --mem=372G
#SBATCH -J hifiasm
#SBATCH -o /staging/biology/u6807805/Tala2/Hifiasm/tala2_hifiasm.o
#SBATCH -e /staging/biology/u6807805/Tala2/Hifiasm/tala2_hifiasm.e

/work/opt/ohpc/Taiwania3/pkg/biology/hifiasm/hifiasm_v0.14/hifiasm \
-o /staging/biology/u6807805/Tala2/Hifiasm/tala2 -t 64 /staging/biology/u6807805/Tala2/m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq
-------------__________________________--------------
> sbatch hifiasm.sh
> squeue -u u6807805
#結果將輸出至/staging/biology/u6807805/Practice/Hifiasm
#後綴bp.p_ctg.gfa為primary assembly 的結果檔


> vim gfa2fasta.sh
-------------_______貼上以下內容________--------------
#!/usr/bin/env bash

if [ -z $1 ]; then
 echo "$0 <gfa file>"
    echo "It will generate the fasta file like ${1%.*}.fasta"
 exit $E_BADARGS
fi
gfatools gfa2fa ${1} > ${1%.*}.fasta
-------------__________________________--------------

#再下指令
> bash gfa2fasta.sh Hifiasm.bp.p_ctg.gfa
#即可在相同資料夾獲得轉換後的fasta files




####Calculate statistics of the assembly####
################# QUAST ################# 
> assembly-stats tala2.bp.p_ctg.fasta
> quast.py tala2.bp.p_ctg.fasta -t 32 -o tala2.bp.p_ctg.fasta.quast
> cd tala2.bp.p_ctg.fasta.quast
> cat report.txt
# quast會生成一個新資料夾，進入資料夾後將report.txt 的內容印在console上，並整理數值至ppt表格中




####### 濾除mitochondiral sequence #######
> wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
## version: 20250905##
> gunzip mitochondrion.1.1.genomic.fna.gz
##跟mitochondiral database進行比對
> minimap2 -x asm5 mitochondrion.1.1.genomic.fna tala2.bp.p_ctg.fasta > mito_all.paf
##禁用 secondary（只保留主比對）
> minimap2 -x asm5 --secondary=no mitochondrion.1.1.genomic.fna tala2.bp.p_ctg.fasta > mito.primary.paf
##限制每條 query 最多回報 1 個 hit
> minimap2 -x asm5 -N 1 mitochondrion.1.1.genomic.fna tala2.bp.p_ctg.fasta > mito.N1.paf

## 濾除50000以下的片段：ptg000009c
> seqkit seq -m 50000 tala2.bp.p_ctg.fasta > tala2.nuclear.fasta
> assembly-stats tala2.nuclear.fasta
> minimap2 -x asm5 mitochondrion.1.1.genomic.fna tala2.nuclear.fasta > mito.paf

#### .paf的看法 ####
## $1  qname     query sequence name        → contig 名稱（ptg000009c）
## $2  qlen      query length               → contig 長度
## $3  qstart    query start (0-based)
## $4  qend      query end
## $5  strand    + or -
## $6  tname     target name                → mitochondrial reference
## $7  tlen      target length
## $8  tstart    target start
## $9  tend      target end
## $10 nmatch    # matching bases
## $11 alnlen    alignment block length
## $12 mapq     mapping quality

> less mito.paf
> cut -f1,2,6,7,10 mito_all.paf | sort -k1,1
##
## $11 = alignment block length
> awk '$11 > 5000' mito_all.paf

##看覆蓋率##
###contigX   contig_length   cov(e.g. 0.95)
> awk '{print $1, $2, $11/$2}' mito_all.paf | sort -k3,3nr
###合併覆蓋率
> awk '{a[$1]+=$11; b[$1]=$2} END {for (i in a) print i, b[i], a[i]/b[i]}' mito_all.paf

##最後用MFannot看ptg000009c是不是粒線體##
https://megasun.bch.umontreal.ca/apps/mfannot/
###詳見mfannot的介紹###




############## QUAST(過濾後) ############## 
> assembly-stats tala2.nuclear.fasta
> quast.py tala2.nuclear.fasta -t 32 -o tala2.nuclear.fasta.quast
> cd tala2.nuclear.fasta.quast
> cat report.txt
# quast會生成一個新資料夾，進入資料夾後將report.txt 的內容印在console上，並整理數值至ppt表格中




######### BUSCO 評估基因體的完整性 #########
##conda activate assembly##
busco -m genome -i [assembly.name].fasta -o [assembly.name] -l ascomycota_odb10 -c 32
> busco -m genome -i tala2.nuclear.fasta -o tala2.nuclear_ascomycota_busco -l ascomycota_odb10 -c 32
> busco -m genome -i tala2.nuclear.fasta -o tala2.nuclear_fungi_busco -l fungi_odb10 -c 32

> busco -m genome -i tala2.nuclear.fasta -o tala2.nuclear_ascomycota_busco_odb12 -l ascomycota_odb12 -c 32
> busco -m genome -i tala2.nuclear.fasta -o tala2.nuclear_fungi_busco_odb12 -l fungi_odb12 -c 32


####用Blobtoolkit篩看看####
#首先建立一個新的environment然後安裝BlobToolKit
#安裝及tutorials 可參考 https://github.com/blobtoolkit/tutorials/tree/main/futurelearn
#相關教學影片 https://www.youtube.com/watch?v=eSUYUQY8r-U&t=2333s&ab_channel=BiodiversityGenomicsAcademy ; https://www.youtube.com/watch?v=SuKBKEH0cMA&t=3500s&ab_channel=BlobToolKit

> conda create -y -n btk python=3.9
> conda activate btk
> pip install blobtoolkit[full]

###Blobtools###
> conda activate btk 
> mkdir /staging/biology/u6807805/Tala2/Blobtools
> mkdir /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir
> cd Blobtools

> mkdir taxdump
#於本地端下載好 https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz 移入上述資料夾，再解壓縮
> tar -xvz -f new_taxdump.tar.gz -C /staging/biology/u6807805/Tala2/Blobtools/taxdump

#進行BlobtoolKit的事前準備
#需準備(1) BUSCO 結果 (2) Read coverage on the assembly (3) Result of BlastN against NCBI nt

####(1) BUSCO analysis to assess the completeness####

##這個前幾步已經完成了###
> conda activate assembly
> cd /staging/biology/u6807805/Tala2//Hifiasm
> busco -m genome -i tala2.bp.p_ctg.fasta -o tala2.bp.p_ctg_ascomycota_busco_odb12 -l ascomycota_odb12 -c 32
> busco -m genome -i tala2.bp.p_ctg.fasta -o tala2.bp.p_ctg_fungi_busco_odb12 -l fungi_odb12 -c 32


#完成後進入所生成之資料夾中，把full_table.tsv (改名以區別basidiomycota 及 fungi的結果) 移到 /staging/biology/u6807805/Practice/Blobtools/HCSR213_blobdir 

-------------------------------
####(2) 計算Read Coverage on the assebled genome####

#將Hifi reads 回貼到assembly上

> conda activate assembly
> minimap2 -ax map-hifi -t 48 /staging/biology/u6807805/Tala2/tala2.nuclear.fasta /staging/biology/u6807805/Tala2/m64531e_251104_190325.bc1001--bc1001.hifi_reads.fastq > tala2_coverage.sam
> samtools sort -@ 32 tala2_coverage.sam -o tala2_coverage.sorted.bam
> samtools depth tala2_coverage.sorted.bam | awk '{sum+=$3} END { print "Average Coverage = ",sum/NR}'

########
samtools faidx tala2.nuclear.fasta 
cut -f1,2 tala2.nuclear.fasta.fai > genome.sizes
bedtools makewindows -g genome.sizes -w 100000 > windows.100kb.bed
bedtools coverage -a windows.100kb.bed -b tala2_coverage.sorted.bam -mean > coverage.100kb.bed
cut -f1-4 coverage.100kb.bed > circos.coverage.100kb.txt 







#########
> samtools index -@ 32 tala2_coverage.sorted.bam tala2_coverage.sorted.bai

#把HCSR213_coverage.sorted.bam移到 /staging/biology/u6807805/Practice/Blobtools/HCSR213_blobdir
-------------------------------
####(3) BlastN against NCBI nt####
#將組裝好的Genome fasta file 使用 blastN blast against NCBI nt database，檢查是否有污染序列
#NCBI nt database 可用電腦至https://ftp.ncbi.nlm.nih.gov/blast/db/ 下載名為 nt.[xxx].tar.gz 之所有檔案，再上傳至國網

> vim blastn.sh
-------------_______貼上以下內容________--------------
#!/bin/bash
#SBATCH -A MST109178
#SBATCH -p ngs53G
#SBATCH -c 8
#SBATCH --mem=53g
#SBATCH -J Blastn
#SBATCH -o /staging/biology/u6807805/Tala2/report/Blastn.o.txt
#SBATCH -e /staging/biology/u6807805/Tala2/report/Blastn.e.txt

/opt/ohpc/Taiwania3/pkg/biology/BLAST/BLAST_v2.15.0/bin/blastn \
-db /staging/biology/u6807805/NCBI_db/nt \
-query /staging/biology/u6807805/Tala2/tala2.nuclear.fasta \
-outfmt "6 qseqid staxids bitscore std" \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-out /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/tala2.bp.p_ctg.out
-------------__________________________--------------
> sbatch blastn.sh
> squeue -u u6807805
> scontrol show job 556537

######看送出的job耗時######
sacct -u u6807805 --starttime 2026-02-09 --format=JobID,JobName,Elapsed,State 

##I. Input file requirements
#Create a metadata file
> vim tala2.yaml
-------------_______貼上以下內容________--------------
assembly:
  alias: tala2
  record_type: contig
taxon:
  name: Talaromyces proteolyticus
-------------__________________________--------------

####Adding input files to blobtoolkit dataset
#adding assembly and a metadata files
> blobtools create --fasta /staging/biology/u6807805/Tala2/tala2.nuclear.fasta --meta tala2.yaml /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/
#會生成 gc.json; identifiers.json; length.json; meta.json; ncount.json 等檔案
#Adding hits
> blobtools add --hits tala2.bp.p_ctg.out --taxrule bestsumorder --taxdump /staging/biology/u6807805/Tala2/Blobtools/taxdump /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/
#會生成許多前綴為bestsumorder的json檔
#Addind coverage 
> blobtools add --cov tala2_coverage.sorted.bam --thread 90 /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/
#Adding BUSCO summary


##到這##
##在xxx_ascomycota_busco_odb12/run_ascomycota_odb12的資料夾裡面 記得要改擋名
> blobtools add --busco tala2_ascomycota_busco_full_table.tsv /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/
> blobtools add --busco tala2_fungi_busco_full_table.tsv /staging/biology/u6807805/Tala2/Blobtools/tala2_blobdir/



###View the dataset on blobtool viewer
#要在電腦local安裝BlobtoolKit，把上述/staging/biology/u6807805/Practice/Blobtools/HCSR213_blobdir/整個資料夾複製到電腦中，再運行以下code，而後將產生的網址複製到任一瀏覽器上
#Teamviewer: 812辦公室電腦(1369561424); 812實驗室電腦(959592747)
#使用812辦公室公用電腦 橘色ubuntu
> conda activate blobtools
> cd /mnt/d/Lucia/Blobtools/
> blobtools view --local HCSR213_blobdir
#注意：若企鵝Ubuntu正在進行運算，會佔用portal，而使blobtools view 無法檢視。
#利用blobtools filter 設定GC threshold 篩選出主基因體
> cd HCSR213_blobdir
> blobtools filter --query-string "http://localhost:8001/view/HCSR213_v2_blobdir/dataset/HCSR213_v2_blobdir/blob?bestsumorder_species--Active=true&catField=bestsumorder_phylum&gc--Min=0.460#Filters" \
--fasta HCSR213.bp.p_ctg.fasta \
/mnt/d/Lucia/Blobtools/HCSR213_blobdir


####Repeatmasker###
# Create the conda environment for RepeatModeler and RepeatMasker
# The RepeatModeler package contains the RepeatMasker package
> conda create -c bioconda -n repeatmodeler
# Activate the environment for RepeatModeler
> conda activate repeatmodeler
# Install the dependencies for RepeatModeler
> conda install -c conda-forge perl-file-which
# Download the NINJA package for large-scale neighbor-joining phylogeny inference and clustering
> cd /staging/biology/u6807805
> mkdir Repeatmasker
> cd Repeatmasker
> wget https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/0.95-cluster_only.tar.gz
> tar -zxvf 0.95-cluster_only.tar.gz
> cd NINJA-0.95-cluster_only/NINJA/
> make # Create the “Ninja” executable file
> pwd
# The pwd Linux command prints the current working directory path
# The standard output of the "pwd" command will be used as a parameter of RepeatModeler
> conda install -c bioconda repeatmodeler -y
> conda install -c bioconda repeatmasker -y


> cd /staging/biology/u6807805/Tala2
> mkdir Repeatmasker
> cd Repeatmasker
> BuildDatabase -name tala2 /staging/biology/u6807805/Tala2/tala2.nuclear.fasta
> RepeatModeler -database tala2 -threads 96 -LTRStruct -ninja_dir /staging/biology/u6807805/Repeatmasker/NINJA-0.95-cluster_only/NINJA  #這步需要較長時間
#-LTRStruct 額外運行LTRharvest和LTR_retriever等工具，尋找LTR element
## HGT ##

###########


> vim repeatmodeler.sh
-------------_______貼上以下內容________--------------
#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J Repeatmodeler_tala22
#SBATCH -p ngs53G
#SBATCH -c 8
#SBATCH --mem=53g
#SBATCH -o RM_tala2.out
#SBATCH -e RM_tala2.err

set -euo pipefail

WORKDIR=/staging/biology/u6807805/repeatmodeler_runs/tala2_${SLURM_JOB_ID}
mkdir -p "$WORKDIR"
cd "$WORKDIR"


source /home/u6807805/miniconda3/etc/profile.d/conda.sh
conda activate repeatmodeler

echo "RepeatModeler:" $(which RepeatModeler)
echo "BuildDatabase:" $(which BuildDatabase)
echo "Perl:" $(which perl)


NINJA_DIR=/staging/biology/u6807805/Repeatmasker/NINJA-0.95-cluster_only/NINJA


GENOME=/staging/biology/u6807805/Tala2/tala2.nuclear.fasta

BuildDatabase -name tala2 "$GENOME"
RepeatModeler -database tala2 -threads ${SLURM_CPUS_PER_TASK} -LTRStruct -ninja_dir "$NINJA_DIR"
-------------__________________________--------------


> sbatch repeatmodeler.sh
> squeue -u u6807805












#NINGA 是一個用於建構大規模親緣關係樹的工具，RepeatModeler在將重複序列分類時需要用到。
> RepeatMasker -lib tala2-families.fa -s -parallel 16 -xsmall -alignments /staging/biology/u6807805/Tala2/tala2.nuclear.fasta -gff
#-s :sensitive mode
#-xsmall :soft masking

#結果會輸出在input 序列的位置，其中 ".tbl" 的檔案紀錄著重複序列統計資訊。".gff" 可作為後續視覺化工具。
> cat tala2.nuclear.fasta.tbl


##Run BRAKER3 on Docker##
##準備gene anoation #將RNA-Seq reads打包,作為BRAKER3預測基因時的參考
##先前的經驗發現CK (對照組)和P（藥劑處理組）的RNA-Seq data預測出來的基因不太一樣，也會影響某些基因的完整度
##最後使用處理組(pencycuron treatment)的RNA-Seq預測基因
> cat 213_P_1_f1_val_1.fq.gz 213_P_2_f1_val_1.fq.gz 213_P_3_f1_val_1.fq.gz > 213_P_cat_f1.fq.gz #把forward reads全部集合成一個檔案 打包
> cat 213_P_1_r2_val_2.fq.gz 213_P_2_r2_val_2.fq.gz 213_P_3_r2_val_2.fq.gz > 213_P_cat_r2.fq.gz #把reverse reads全部集合成一個檔案
> gunzip 213_P_cat_f1.fq.gz
> gunzip 213_P_cat_r2.fq.gz
##將 213_P_cat_f1.fq.gz 和 213_P_cat_r2.fq.gz 下載至812辦公室之公用電腦 /mnt/d/Lucia/BRAKER_gene_prediction/HCSR213/HCSR213RNAseq

#下載OrthoDB之Fungi.fa.gz (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)

#用812公用電腦執行，打開企鵝圖標的ubuntu 以及Docker等應用程式
##Mount data path onto the container
> docker run --user 1000:100 --rm -it -v /mnt/d/Lucia/BRAKER_gene_prediction/HCSR213_v2/:/home/jovyan teambraker/braker3:latest bash
##Commands for BRAKER3
> braker.pl --genome=HCSR213.bp.p_ctg.filtered.fasta.masked.fasta \
       --rnaseq_sets_ids=HCSR213_PNC_cat_f1,HCSR213_PNC_cat_r2 \
       --rnaseq_sets_dirs=HCSR213RNAseq \
       --prot_seq=Fungi.fa

#--genome=將進行基因預測之基因體檔案名稱
#--rnaseq_sets_ids=RNA-Seq data之檔案名稱(不加副檔名)
#--rnaseq_sets_dirs=儲存RNA-Seq data之相對路徑


####無RNA-Seq data之案例
#利用近緣種之proteins來輔助基因預測
#從NCBI 下載近緣種基因體之protein.fasta，亦可與OrthoDB的Fungi.fa打包成一個檔案，用於輔助預測
> docker run --user 1000:100 --rm -it -v /mnt/d/KCC/BRAKER_gene_prediction/tala2/:/home/jovyan teambraker/braker3:latest bash
> braker.pl \
  --genome=tala2.nuclear.fasta \
  --prot_seq=Fungi.fa \
  --AUGUSTUS_ab_initio


#送出指令後會顯示 # Fri Dec  5 09:14:10 2025: Log information is stored in file /home/jovyan/braker/braker.log 的字樣，並閃爍游標。屬正常現象。
#打開資料夾 HCSR213 會發現多了一個叫braker的資料夾
#點進去可以打開 braker.log這個檔案，會顯示當前分析的進度。想要更新內容的話，須先關閉再重新打開檔案。
##大約需要16-25小時的運行時間，取決於基因體大小，運行後會獲得 (1)xxx.gtf (2)xxx.codingseq (3)xxx.aa
#將xxx.aa 用於後續基因功能註釋分析

###BUSCO 檢測預測出來的基因的完整度###
> busco --list-datasets
> busco -m proteins -i HCSR213.braker.aa -o HCSR213.braker.aa_busco_basidiomycota_odb10 -l basidiomycota_odb10 -c 32

busco -m proteins -i tala2.braker.aa -o tala.braker.aa_busco_ascomycota_odb12 -l ascomycota_odb12 -c 32
busco -m proteins -i tala2.braker.aa -o tala.braker.aa_busco_eurotiomycetes_odb12 -l eurotiomycetes_odb12 -c 32
busco -m proteins -i tala2.braker.aa -o tala.braker.aa_busco_eurotiales_odb12 -l eurotiales_odb12 -c 32
busco -m proteins -i tala2.braker.aa -o tala.braker.aa_busco_aspergillus_odb12 -l aspergillus_odb12 -c 32
busco -m proteins -i tala2.braker.aa -o tala.braker.aa_busco_penicillium_odb12 -l penicillium_odb12 -c 32



-------------------------------------------------------------------------------------------------------
####Gene Functional Annoatation####
#參考教學：https://yanzhongsino.github.io/2021/05/17/omics_genome.annotation_function/

###1. NCBI nr
###把xxx.aa拿去blastp NCBI nr database###
###Database下載：https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz (~186G)，下載後上傳到國網用Dimand 製作成他適用的的database
#Diamond make NR database#
#之前已經建過了，可以不用再一次
> vim diamond_makedb.sh
-------------_______貼上以下內容________--------------
#!/bin/bash
#!/bin/bash
#SBATCH -A MST109178
#SBATCH -p ngs53G
#SBATCH -c 8
#SBATCH --mem=53g
#SBATCH -J Diamond
#SBATCH -o /staging/biology/u6807805/Tala2/report/Diamond.mkdb.o.txt
#SBATCH -e /staging/biology/u6807805/Tala2/report/Diamond.mkdb.e.txt

/opt/ohpc/Taiwania3/pkg/biology/diamond/diamond_v2.0.11/diamond \
makedb --in /staging/biology/u6807805/nr --db nr
-------------__________________________--------------

#再將xxx.aa.fasta 拿去blast nr database 進行註釋#
> vim blastnr.sh
-------------_______貼上以下內容________--------------
#!/bin/bash
#SBATCH -A MST109178
#SBATCH -p ngs53G
#SBATCH -c 8
#SBATCH --mem=53g
#SBATCH -J Diamond
#SBATCH -o /staging/biology/u6807805/Tala2/report/Diamond.blastp.o.txt
#SBATCH -e /staging/biology/u6807805/Tala2/report/Diamond.blastp.e.txt
/opt/ohpc/Taiwania3/pkg/biology/diamond/diamond_v2.0.11/diamond \
blastp --db /staging/biology/u6807805/nr.dmnd --query /staging/biology/u6807805/Tala2/Braker/tala2.braker.aa.fasta \
-o /staging/biology/u6807805/Tala2/tala2.nr.protein_annoation.txt --evalue 0.00001 \
--max-target-seqs 5 --outfmt 6
-------------__________________________--------------
> sbatch blastnr.sh
> squeue -u u6807805
#evalue




###2. InterProScan5
##!!!!!要記得先去除xxx.aa.fasta 中每個amino acid sequence末端的* 不然會馬上報錯！可使用Subline 將其打開，用Replace 功能。
##InterProScan5安裝於812公用電腦d槽，進入該資料夾後 /mnt/d/Interproscan/interproscan-5.68-100.0
#執行
> ./interproscan.sh -i /mnt/d/Lucia/BRAKER_gene_prediction/HCSR213/braker/HCSR213_braker.aa.fasta -b HCSR213_PNC_.out.iprscan -goterms -iprlookup -pa -dp -cpu 24

#會看到以下進度，直到它顯示InterProScan analyses completed
15/07/2024 11:11:01:436 Welcome to InterProScan-5.68-100.0
15/07/2024 11:11:01:436 Running InterProScan v5 in STANDALONE mode... on Linux
15/07/2024 11:11:08:022 RunID: DESKTOP-7J5L5BJ_20240715_111106969_h15x
15/07/2024 11:11:15:236 Loading file /mnt/d/HCSR214.final.braker.aa.04.fasta
15/07/2024 11:11:15:243 Running the following analyses:
[AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_05,MobiDBLite-2.0,NCBIfam-14.0,PANTHER-18.0,Pfam-37.0,PIRSF-3.10,PIRSR-2023_05,PRINTS-42.0,ProSitePatterns-2023_05,ProSiteProfiles-2023_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75]
Pre-calculated match lookup service DISABLED.  Please wait for match calculations to complete...
15/07/2024 12:38:45:450 25% completed
15/07/2024 13:11:28:201 50% completed
15/07/2024 13:43:42:866 75% completed
15/07/2024 14:43:42:878 87% completed
15/07/2024 15:43:42:977 91% completed
15/07/2024 16:43:43:035 91% completed
15/07/2024 17:43:43:155 91% completed
15/07/2024 18:43:43:258 92% completed
15/07/2024 19:43:43:351 92% completed
15/07/2024 20:43:43:435 92% completed
15/07/2024 21:43:43:534 98% completed
15/07/2024 22:43:43:633 98% completed
15/07/2024 23:30:25:437 100% done:  InterProScan analyses completed

##完成後取 *.out.tsv檔來進行後續資料整理
##註：由於HCSR214蛋白數太多，電腦會crush，因此將其拆分成四個檔案分別執行，最後再把tsv檔案整合


###3. eggNOG-mapper
#可注釋GO, KEGG等
#webtool: http://eggnog-mapper.embl.de/
#上傳xxx.aa.fasta使用預設參數即可
#完成後會有email通知，下載out.emapper.annotations，改檔名為.txt，再做後續資料整理

###4. Pannzer2
##此用Uniport database 進行註釋，可預測GO 及產出Functional description
#webtool: http://ekhidna2.biocenter.helsinki.fi/sanspanz/
#上傳xxx.aa.fasta使用預設參數即可
#完成後會有email通知，結果*.anno.out.txt; *.GO.out.txt 打開頁面後，右鍵儲存
http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//vNEDvStZnt1/index.html


###5. KAAS
#可註釋KEGG
#webtool: https://www.genome.jp/kaas-bin/kaas_main
#上傳xxx.aa.fasta，並選擇GENES data set
cne, cgi, tms, ppl, tvs, hir, mpr, scm, uma, mgl, mrt, sce, ago, kla, kmx, ncs, tpf, pic, cal, pkz, yli, ncr, mgr, fgr, fpu, fox, nhe, maw, ssl, bfu, psco, bze, pno, ani, afm, aor, afv, pcs, ang
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1766381214&key=PItF9KYv




------------------------------------
#####將上述資料進行整理，並合併#####
#從interproscan結果提取interpro的數據
> cat *.iprscan.tsv |cut -f1,12,13|grep "IPR"|awk '{print $1"\t"$2": "$3}'|sort -k 1.1n|uniq |awk -F"\t" '{a[$1]=a[$1]$2"; "}END{for( i in a){print i"\t"a[i]}}'|sed "s/; $//g" |sort -k 1.1n |sed '1i\gene\tinterpro' >tala2_PNC_iprs.interpro.anno
#從interproscan結果提取pfam的數據
> cat *.iprscan.tsv |awk '$4 == "Pfam" {print $1"\t"$5": "$6}' |sort -k 1.1n |uniq|awk -F"\t" '{a[$1]=a[$1]$2"; "}END{for( i in a){print i"\t"a[i]}}'|sed "s/; $//g" |sort -k 1.1n |sed '1i\gene\tinterproscan.pfam'>tala2_PNC_iprs.pfam.anno
#從interproscan結果提取PANTHER的數據
> cat *.iprscan.tsv |awk '$4 == "PANTHER" {print $1"\t"$5": "$6}'|sort -k 1.1n|uniq |awk -F"\t" '{a[$1]=a[$1]$2"; "}END{for( i in a){print i"\t"a[i]}}'|sed "s/; $//g" |sort -k 1.1n |sed '1i\gene\tinterproscan.panther'>tala2_PNC_iprs.panther.anno


#從Pannzer2結果提取uniport的數據
> cat tala2.anno.out.txt |awk '$2 == "DE" {print $1"\t"$6}' |sort -k 1.1n |uniq |sed '1i\gene\tpannzer2.uniprot' >tala2_PNC_pannzer2.uniprot.anno
#從Pannzer2結果提取GO的數據
#出來的data在用,取代"; "
> cat tala2.GO.out.txt |awk '{print $1"\tGO:"$3}' |awk -F"\t" '{a[$1]=a[$1]$2"; "}END{for( i in a){print i"\t"a[i]}}'|sed "s/; $//g" |sort -k 1.1n |uniq |sed "1s/.*/gene\tpannzer2.go/" >tala2_PNC_pannzer2.go.anno

#從eggNOG-mapper結果提取數據
> cat tala2.out.emapper.annotations |sed '/#.*/d' |cut -f1,8,9 |awk -F"\t" '{print $1"\t"$3": "$2}' |sed '/-: -/d'|sort -k 1.1n |uniq |sed '1i\gene\teggnog.preferred.name: description' >tala2_PNC_eggnog.anno
> cat tala2.out.emapper.annotations |sed '/#.*/d' |cut -f1,10 |sed '/\t-$/d' |sort -k 1.1n |uniq |sed '1i\gene\teggnog.go' >tala2_PNC_eggnog.go.anno
> cat tala2.out.emapper.annotations |sed '/#.*/d' |cut -f1,12 |sed '/\t-$/d' |sort -k 1.1n |uniq |sed '1i\gene\teggnog.kegg' >tala2_PNC_eggnog.kegg.anno

#將結果檔案轉成.txt 並修改geneID標題為GeneID
#匯入R進行合併 
#首先利用aa.fasta 提取基因id (Genome_assembly/Final_assembly/Extract Gene ID from fasta.R)，並用merge function把其餘的結果整合上去，最後將GO ID整合 (Genome_assembly/Final_assembly/213Merge_data.R)
#最後匯出為 (1)Genome_assembly/Final_assembly/tala2_final/Functional_annoation/tala2_PNC_merged_functional_annotations.txt (2)Genome_assembly/Final_assembly/HCSR214_final/Functional_annoation/HCSR214_PNC_merged_functional_annotations.txt


############################
####完成基因體組裝及基因註釋####
############################






Ctrl + Z
>>>[1]+  Stopped   python main.py


bg
>>>[1]+ python main.py &


jobs -l
disown -h %1

ps -u $USER | grep python
ps -u $USER | grep Rscript
ps -fp <PID>






#######找transposon#######
> conda activate repeatmodeler
> cd /staging/biology/u6807805/Tala2/transposon
> BuildDatabase -name tala2_trans /staging/biology/u6807805/Tala2/transposon/flanking433.fasta
> RepeatModeler -database tala2_trans -threads 96 -LTRStruct -ninja_dir /staging/biology/u6807805/Repeatmasker/NINJA-0.95-cluster_only/NINJA 

#NINGA 是一個用於建構大規模親緣關係樹的工具，RepeatModeler在將重複序列分類時需要用到。
> RepeatMasker -species Fungi -s -parallel 16 -xsmall -alignments /staging/biology/u6807805/Tala2/transposon/flanking433.fasta -gff

> RepeatMasker -lib tala2-families.fa -s -parallel 16 -xsmall -alignments /staging/biology/u6807805/Tala2/transposon/flanking433.fasta -gff
#-s :sensitive mode
#-xsmall :soft masking

#結果會輸出在input 序列的位置，其中 ".tbl" 的檔案紀錄著重複序列統計資訊。".gff" 可作為後續視覺化工具。
> cat tala2.nuclear.fasta.tbl
