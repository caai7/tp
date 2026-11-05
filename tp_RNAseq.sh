### Tp RNA-seq ###

#---1-- 檢查Read品質 ---#
# Raw FASTQ -> FastQC (檢查read品質) -> MultiQC (整合QC報告) -> fastp (依照結果決定 移除adapter或修剪低品質reads) 
# -> FastQC/MultiQC 重新檢查品質 -> Alignment

##檢查檔案大小
ls -lh *.fastq.gz




#---2-- 比對到 Reference genome ---#
