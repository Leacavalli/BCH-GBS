
# Code
### 0. Navigate to file location

cd /n/netscratch/hanage_lab/Lab/lcavalli/BCH_GBS/Pyseer

### 1. Install Pyseer Files

cd Files
git clone https://github.com/mgalardini/pyseer

### 2. Install Pyseer through conda
module load python
# conda create -n pyseer
conda activate pyseer
# conda install pyseer
pyseer -h

### 3. Run Pyseer for each phenotype of interest
# cd /n/netscratch/hanage_lab/Lab/lcavalli/BCH_GBS/Pyseer/ICU_Infants
sbatch -p shared -c 8 -t 0-00:10 --mem=10000 --wrap="python /n/netscratch/hanage_lab/Lab/lcavalli/BCH_GBS/Pyseer/pyseer/scripts/phylogeny_distance.py --lmm tree_GBS_BCH.tre > phylogeny_K.tsv"
sbatch -p test -c 8 -t 0-00:10 --mem=10000 --wrap="pyseer --lmm --phenotypes phenotypes.txt --pres gene_presence_absence.txt --similarity phylogeny_K.tsv --cpu 8 > Output_GAS.txt"


### 3.1. Run Pyseer for continuous age
# cd /n/netscratch/hanage_lab/Lab/lcavalli/BCH_GBS/Pyseer/LOD_vs_VLOD
sbatch -p test -c 8 -t 0-00:10 --mem=10000 --wrap="pyseer --lmm --phenotypes continuous_phenotypes.txt --pres gene_presence_absence.txt --similarity phylogeny_K.tsv --cpu 8 > Output_GAS_continuous.txt"
