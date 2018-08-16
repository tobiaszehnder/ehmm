#! /bin/bash

# make sure a launcher called `ehmm.R` has been created using the `ehmm:::getLauncher` function
# make sure the launcher is in a folder listed in your `PATH` environment variable (`export PATH=/path/to/the/launcher/:$PATH`)

# pass number of threads (ideally 21, one for each region)
if [ $# -ne 1 ]; then echo "Usage: ./script.sh <nthreads>"; exit 1; fi 

nthreads=$1
root=$(realpath .)
indir=$root/inputfiles/
bamdir=$indir/bam
outdir=$root/output/

mkdir -p $bamdir $outdir

# download and preprocess bam files
make

# download genome file and create regions file
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O $indir/mm10.genome
cat $indir/mm10.genome | grep -v "M" | grep -v "Un" | grep -v "random" | awk 'BEGIN{FS=OFS="\t"} {print $1, 100, int($2/100)*100}' > $indir/regions.bed 

# run eHMM
echo 'run eHMM'
ehmm.R applyModel \
       --regions $indir/regions.bed \
       --genomeSize $indir/mm10.genome \
       --mark ATAC-seq:$bamdir/ATAC-seq.bam \
       --mark H3K27AC:$bamdir/H3K27AC.bam \
       --mark H3K4ME1:$bamdir/H3K4ME1.bam \
       --mark H3K4ME3:$bamdir/H3K4ME3.bam \
       --provideModel \
       --nthreads $nthreads \
       --outdir $outdir
