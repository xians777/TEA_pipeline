# Pipeline for calling TE insertion

##Dependency (configuration on O2)
+ module load gcc/6.2.0
+ module load boost/1.62.0
+ module load pigz/2.3.4
+ module load R/3.4.1
+ module load samtools/1.3.1
+ module load bedtools/2.26.0

##Run the pipeline

###Run the calling step
```
sh 00_main-1filelist-2title-3cores-4mem.sh input_bams.list output_folder 8 100
```
###Run the filtering step

