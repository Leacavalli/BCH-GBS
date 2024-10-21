

* Make sure that your raw reads are in 0.RAW_READS/ and called S*_{1,2}.fastq.gz
* Make sure that the Outgroup for the phylogeny is in 0.RAW_READS/ and called Outgroup_{1,2}.fastq.gz

# 1. Run Trimmomatic
### 1.1. Install Software and dependencies
```
module load python
# conda create -n java
source activate java
# conda install openjdk=11
```
### 1.2. Install dependencies
```
mkdir Files
cd Files
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar -h
```
### 1.3. Run
```
mkdir 1.TRIMMED_READS
cd 1.TRIMMED_READS
for i in ../0.RAW_READS/*1.fastq.gz
do
 base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="java -jar ../Files/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog S1_1.log ../0.RAW_READS/$base\_1.fastq.gz ../0.RAW_READS/$base\_2.fastq.gz $base\_trimmed_1.fq.gz $base\_unpaired_1.fq.gz $base\_trimmed_2.fq.gz $base\_unpaired_2.fq.gz ILLUMINACLIP:../Files/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50 > $base\_trim.out"
done
cd ../
```

# X. (OPTIONAL) Run Kraken on trimmed reads
# This is optional because FastANI is later ran on assemmblies
### 2.1. Run
```
mkdir 2.KRAKEN
cd 2.KRAKEN
for i in ../1.TRIMMED_READS/*trimmed*.fq.gz
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
  sbatch -p shared -t 0-00:10 --mem=100000 --wrap="kraken2 --db /n/holylfs05/LABS/hanage_lab/Lab/holyscratch01/lcavalli/Software/kraken2/standard_db --report $base\_report.txt --paired ../1.TRIMMED_READS/$base\_trimmed_1.fq.gz ../1.TRIMMED_READS/$base\_trimmed_2.fq.gz > $base\_kraken.out" ;
done
```
### 2.2. Summarise
```
echo -e "Sample\tPercent GBS Reads" > KRAKEN_contaminated.txt
for i in ../1.TRIMMED_READS/*trimmed*.fq.gz
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
  percent=$(grep 'Streptococcus agalactiae' $base\_report.txt | head -1 | awk '{print $1}')
  echo -e $base'\t'$percent
done  >> KRAKEN_contaminated.txt
```
# 2. Run FASTQC on trimmed reads
### 2.1. Install Software and dependencies
```
module load python
# conda create -n fastqc
source activate fastqc
# conda install fastqc
```
### 2.2. Run
```
mkdir 2.FASTQC
cd 2.FASTQC
for i in ../1.TRIMMED_READS/*trimmed*.fq.gz
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-3)
  sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="fastqc $i --extract -o ."
done
# Check results
cat */summary.txt | grep 'FAIL'
cd ../
```

# 3. Run srst2 on trimmed reads
### 3.1. Install Software and dependencies
```
module load python
# conda create -n srst2
source activate srst2
# conda install srst2
```
### 3.2. Run
```
mkdir 3.SEROTYPE_MLST
cd 3.SEROTYPE_MLST
for i in ../1.TRIMMED_READS/*_trimmed_1.fq.gz
do
 base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
 sbatch -p shared -c 1 -t 0-00:30 --mem=10000 --wrap="srst2 --input_pe ../1.TRIMMED_READS/$base\_trimmed_{1,2}.fq.gz --output $base --log --gene_db ../Files/GBS-SBG.fasta --mlst_db ../Files/Streptococcus_agalactiae.fasta --mlst_definitions ../Files/profiles_csv --mlst_delimiter '_'"
done
cd ../
```
### 3.#. Summarise Results
```
# SEROTYPE
cd 3.SEROTYPE_MLST/
awk 'FNR==1 && NR!=1 {next} {print}' *__genes__GBS-SBG__results.txt > BCH_GBS_serotype.txt
# ST
cd 3.SEROTYPE_MLST/
awk 'FNR==1 && NR!=1 {next} {print}' *_mlst__Streptococcus_agalactiae__results.txt > BCH_GBS_mlst.txt
```

# 4. Run Sanger typer
## Note: This needs to be run on my laptop, rather than the RC because it requires Docker.
## 4.0. Clone repository
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```
## 4.1. Run only Surface typer
```
for i in *_trimmed_1.fq.gz
do
 base=$(echo $i | cut -d '/' -f8-8 | cut -d '_' -f1-1)
 nextflow run main.nf --reads $base\_trimmed_{1,2}.fq.gz --results_dir $base\_surfacetyper_results --run_sero_res false --run_mlst false
done
```
## 4.2. Run only serotyping and resistance workflows
```
for i in *_trimmed_1.fq.gz
do
 base=$(echo $i | cut -d '/' -f8-8 | cut -d '_' -f1-1)
 nextflow run main.nf --reads $base\_trimmed_{1,2}.fq.gz --results_dir $base\_sero_res_results --run_surfacetyper false --run_mlst false
done
```
# 4.3. Join results
```
# SURFACE PROTEINS
awk 'FNR==1 && NR!=1 {next} {print}' *_surfacetyper_results/surface_protein_incidence.txt> BCH_GBS_surface_protein.txt
# RESISTANCE MUTATIONS
awk 'FNR==1 && NR!=1 {next} {print}' *_sero_res_results/gbs_res_variants.txt > BCH_GBS_Resistance_MUTS.txt
# RESISTANCE GENES
awk 'FNR==1 && NR!=1 {next} {print}' *_sero_res_results/serotype_res_incidence.txt > BCH_GBS_Resistance_Genes.txt
```

# 4. Run SNIPPY
### 4.1. Install Software and dependencies
```
module load python
# conda create -n snippy
source activate snippy
# conda install snippy
```
### 4.2. Run
```
mkdir 4.SNIPPY/
cd 4.SNIPPY/
for i in ../1.TRIMMED_READS/*_trimmed_1.fq.gz
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="snippy --outdir ../4.SNIPPY/$base --R1 ../1.TRIMMED_READS/$base\_trimmed_1.fq.gz --R2 ../1.TRIMMED_READS/$base\_trimmed_2.fq.gz --ref ../Files/AP018935.1.fa"
done
cd ../
```


# 5. Run ASSEMBLY
### 5.1. Install Software and dependencies
```
module load python
# conda create -n unicycler
source activate unicycler
# conda install bioconda::unicycler
```
### 5.2. Run
```
mkdir 5.ASSEMBLY/
cd 5.ASSEMBLY/
for i in ../1.TRIMMED_READS/*_trimmed_1.fq.gz
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '_' -f1-1)
 sbatch -p shared -c 1 -t 0-02:00 --mem=10000 --wrap="unicycler -1 ../1.TRIMMED_READS/$base\_trimmed_1.fq.gz  -2 ../1.TRIMMED_READS/$base\_trimmed_2.fq.gz -o $base\.unicycler.out
 mv $base\.unicycler.out/assembly.fasta $base\.fasta"
done
cd ../
```

# 6. Run QUAST
### 6.1. Install Software and dependencies
```
module load python
# conda create -n quast
source activate quast
# conda install quast
```
### 6.2. Run
```
mkdir 6.QUAST
cd 6.QUAST
for i in ../5.ASSEMBLY/*fasta
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '.' -f1-1)
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="quast -o $base\.quast.out $i"
done
cd ../
```

# 7. Run FastANI
### 7.1. Install Software and dependencies
```
module load python
# conda create -n fastani
source activate fastani
# conda install fastani
```
### 7.2. Run
```
mkdir 7.FastANI
cd 7.FastANI
for i in ../5.ASSEMBLY/*fasta
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '.' -f1-1)
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="fastANI -q $i -r ../Files/AP018935.1.fa -o $base\.fastani.out"
done
cd ../
```


# 8. Run ABRICATE
### 8.1. Install Software and dependencies
```
module load python
# conda create -n abricate
source activate abricate
# conda install abricate
```
### 8.2. Run
```
mkdir 8.ABRICATE
cd 8.ABRICATE
for i in ../5.ASSEMBLY/*fasta
do
  base=$(echo $i | cut -d '/' -f3-3 | cut -d '.' -f1-1)
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="abricate --db vfdb $i > $base\_vfdb.tab"
 sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="abricate --db sip_gene $i > $base\_sip_gene.tab"
done
cd ../
```

# 0. QC STEPS: Filter out low quality isolates
```
## Create necessary directories if they don't already exist
  mkdir -p 4.SNIPPY/Filter_out_QC
  mkdir -p 5.ASSEMBLIES/Filter_out_QC
  mkdir -p 8.ABRICATE/Filter_out_QC

  for i in 0.RAW_READS/*1.fastq.gz
do
  base=$(echo $i | cut -d '/' -f2 | cut -d '_' -f1)

  # Extract information from the QUAST report
  ANI=$(awk '{print $3}' 7.FastANI/$base.fastani.out)
  n50=$(grep -w "N50" 6.QUAST/$base.quast.out/report.txt | awk '{print $2}')
  gc=$(grep -w "GC (%)" 6.QUAST/$base.quast.out/report.txt | awk '{print $3}')
  total_length=$(grep -w "Total length" 6.QUAST/$base.quast.out/report.txt | grep -v "(>= " | awk '{print $3}')
  contigs=$(grep "^# contigs " 6.QUAST/$base.quast.out/report.txt | tail -1 | awk '{print $3}')
  Adapter_content_1=$(grep "Adapter Content" 2.FASTQC/${base}_trimmed_1_fastqc/summary.txt | awk '{print $1}')
  Adapter_content_2=$(grep "Adapter Content" 2.FASTQC/${base}_trimmed_2_fastqc/summary.txt | awk '{print $1}')

  # Check the conditions and move files if necessary
  if (( $(echo "$ANI < 95" | bc -l) )) ||
     [ "$n50" -lt 30000 ] ||
     [ "$contigs" -gt 250 ] ||
     (( $(echo "$gc <= 30" | bc -l) )) ||
     (( $(echo "$gc >= 40" | bc -l) )) ||
     [ "$total_length" -lt 1700000 ] ||
     [ "$total_length" -gt 2400000 ] ||
     [[ "$Adapter_content_1" == "FAIL" || "$Adapter_content_2" == "FAIL" ]]
  then
    # Check if the SNIPPY directory exists, and move it if it does
    if [ -d "4.SNIPPY/$base" ]; then
      mv 4.SNIPPY/$base 4.SNIPPY/Filter_out_QC/
    fi

    # Check if the assembly fasta file exists, and move it if it does
    if [ -f "5.ASSEMBLIES/$base.fasta" ]; then
      mv 5.ASSEMBLIES/$base.fasta 5.ASSEMBLIES/Filter_out_QC/
    fi

    # Check if the abricate output file exists, and move it if it does
    if [ -f "8.ABRICATE/${base}_vfdb.tab" ]; then
      mv 8.ABRICATE/${base}_vfdb.tab 8.ABRICATE/Filter_out_QC/
    fi

    if [ -f "8.ABRICATE/${base}_sip_gene.tab" ]; then
      mv 8.ABRICATE/${base}_sip_gene.tab 8.ABRICATE/Filter_out_QC/
    fi
  fi
done

```



# 9. ASSIGN_CC
```
# Create Output Directory
mkdir 9.ASSIGN_CC

# Append all results from mlst to a single file
awk '(NR == 1) || (FNR > 1)' 3.SEROTYPE_MLST/*__mlst__Streptococcus_agalactiae__results.txt > 3.SEROTYPE_MLST/combined_results.txt
cut -f1,2 3.SEROTYPE_MLST/combined_results.txt  | tail -n +2  | sed 's/\t/,/g' | sort -t',' -k2,2 | sed 's/\*//g' > 9.ASSIGN_CC/STs.txt

# Define the file paths
cut -d',' -f1,9 Files/GBS_CC_profiles.csv | tail -n +2  | sort -t',' -k1,1 > 9.ASSIGN_CC/GBS_ST_CC_ref.txt

# Process the files and join them based on the 'ST' column
echo "ST,Sample,clonal_complex" > 9.ASSIGN_CC/combined_results_with_cc.txt
join -t',' -1 2 -2 1 -a 1 -e 'NA' -o 1.2,1.1,2.2 9.ASSIGN_CC/STs.txt 9.ASSIGN_CC/GBS_ST_CC_ref.txt -a1  >> 9.ASSIGN_CC/combined_results_with_cc.txt
```

# 10. ABRICATE_SUMMARY
### 10.1. Install Software and dependencies
```
module load python
# conda create -n abricate
source activate abricate
# conda install abricate
```
### 10.2. Summarise results
```
cd 8.ABRICATE
sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="abricate --summary *_vfdb.tab > Summary_vfdb.tab"
sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="abricate --summary *_sip_gene.tab > Summary_sip_gene.tab"
```

# 11. SNIPPY_MULTI
### 11.1. Install Software and dependencies
```
module load python
# conda create -n snippy
source activate snippy
# conda install snippy
```
### 11.2. Summarise results
```
mkdir 10.SNIPPY_MULTI
cd 10.SNIPPY_MULTI
sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="snippy-core --ref ../Files/AP018935.1.fa --prefix core_outgroup ../4.SNIPPY/S* ../4.SNIPPY/Outgroup"
sbatch -p shared -c 1 -t 0-00:10 --mem=10000 --wrap="snippy-clean_full_aln core_outgroup.full.aln > clean_outgroup.full.aln"
```


# 12. RAxML
### 12.1. Install Software and dependencies
```
cd Files
wget https://github.com/stamatak/standard-RAxML/archive/master.zip
unzip master.zip
# compile source code
cd standard-RAxML-master/
make -f Makefile.SSE3.PTHREADS.gcc # parallelized and x86 processor optimized version
cd ../
standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -h
```
### 12.2. Summarise results
```
mkdir 12.RAxML
cd 12.RAxML
# Generate 100 ML trees on distinct starting trees and output the best likelihood tree
sbatch -p test -c 50 -t 0-12:00 --mem=100000 --wrap="/n/holyscratch01/hanage_lab/lcavalli/BCH_GBS/NEXTFLOW_PIPELINE/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -# 100 -s ../11.SNIPPY_MULTI/clean_outgroup.full.aln -n T1 -o Outgroup"
# Generate 250 bootstrap tree to infer statistical support of the branches.
sbatch -p test -c 50 -t 0-12:00 --mem=100000 --wrap="/n/holyscratch01/hanage_lab/lcavalli/BCH_GBS/NEXTFLOW_PIPELINE/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRGAMMA -p 12345 -b 12345 -# 250 -s ../11.SNIPPY_MULTI/clean_outgroup.full.aln -n T2 -o Outgroup"
# Draw bipartitions on the best ML tree
sbatch -p test -c 50 -t 0-12:00 --mem=100000 --wrap="/n/holyscratch01/hanage_lab/lcavalli/BCH_GBS/NEXTFLOW_PIPELINE/Files/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 50 -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T1 -z RAxML_bootstrap.T2 -n T3"
```



#  MASHTREE
```
module load python
source activate mashtree
sbatch -p shared -t 0-01:00 --mem=10000 --wrap="mashtree *.fastq* > tree3.dnd"
```
