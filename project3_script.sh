#!/bin/bash

#Please use absolute path for input_directory
input_directory="$1"

### ENVIRONMENT CREATION

#creation of MLST environment
conda create -n mlst -c conda-forge -c bioconda mlst -y

#creation of maast environment and activation of environment
conda create -n maast python=3.12
conda activate maast
#Required libraries and programs
conda install -c conda-forge numpy=1.19.5
conda install -c conda-forge scipy=1.5.4
conda install -c conda-forge biopython=1.83
conda install -c conda-forge networkx=2.5.1
conda install -c conda-forge scipy=1.12.0
conda install -c bioconda mash mummer
conda install -c bioconda mummer4
conda install -c bioconda pigz lbzip2 lz4
conda install -c bioconda fasttree
conda deactivate maast

#creation of checkm environment
conda create -n checkm -c conda-forge -c bioconda checkm-genome -y

# LisSero

conda create --name lissero
conda activate lissero
pip3 install lissero
conda install -c bioconda blast=2.10.1
conda deactivate lissero

# Quast
conda create -n qual_eval python=3.7.12 -y
pip install quast

### RUN COMMANDS

#MLST COMMAND

conda activate mlst
mkdir ./mlst
#copies raw data into form suitable for MLST
cp -r "$input_directory"/* ./mlst/
cd mlst
#Run mlst
mlst *.fna > MLST_Summary.tsv
cd ..
conda deactivate


#MAAST COMMANDS
conda activate maast
#Copied a copy of Maast to my local computer
git clone https://github.com/zjshi/Maast.git
cd Maast # Navigating to directory where maast is present
make #This compiles the source code of maast
chmod 755 maast #to make GT-Pro ready to execute
#Maast command for genotyping

mkdir ./maast_output
maast end_to_end --in-dir $input_directory --out-dir ./maast_output --min-prev 0.9 --snp-freq 0.01
#the step generates a list of input pairs. 
paste <(find ./maast_output/gt_results/ -name '*.tsv' | sort) <(find ./maast_output/gt_results/ -name '*.tsv' | sort | cut -d'/' -f4 | cut -d'.' -f1) > genotypes.input.tsv
./Maast tree --input-list ./genotypes.input.tsv --out-dir ./tree_results/ #Phylogenetic tree building
cd ..
conda deactivate


#KRAKEN2 COMMANDS
# install kraken2
git clone https://github.com/DerrickWood/kraken2.git
./install_kraken2.sh $KRAKEN2_DIR # $KRAKEN2_DIR: directory where you want to install Kraken
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin
 
# install DB for kraken2 (need ~100GB)
kraken2-build --standard --threads 24 --db $DBNAME # $DBNAME: your preferred location
 
files=$(find "$input_directory" -type f -name "*.fa")
 
for file in "$input_directory"/*.fa; do
  sample_id=$(basename "$file" .fa)
  output_path="./kraken-out/${sample_id}.out"
  report_path="./kraken-report-short/${sample_id}.report"
  kraken2 --db Standard-8 $file --threads 128 --output $output_path --report $report_path
done


#LISSERO

conda activate lissero
mkdir lissero_output

# Create an empty file to store the concatenated results
touch lissero_output/all_lissero_results.txt

# Iterate over each .fa file in the filtered_contigs directory
for fasta_file in "$input_directory"/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}" 

    # Run lissero on the current fasta file and store the output in a temporary file
    lissero "$fasta_file" > "lissero_output/${filename_no_ext}_lissero_results.txt.tmp"
    
    # Append the contents of the temporary file to the concatenated results file
    cat "lissero_output/${filename_no_ext}_lissero_results.txt.tmp" >> lissero_output/all_lissero_results.txt
    
    # Remove the temporary file
    rm "lissero_output/${filename_no_ext}_lissero_results.txt.tmp"
done
conda deactivate
#BLASTN

cd filtered_contigs_by_serotype
mkdir blast_out_serotype1_2a
### Run blastn on all the fasta files and output in format 6
### Use respective reference genomes of each classified serotype from LisSero
for fasta_file in serotype1_2a/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}"

    blastn -query "$fasta_file" -subject GCF_000196035.1_ASM19603v1_genomic.fna \
        -out "blast_out_serotype1_2a/${filename_no_ext}_serotype1_2a_blast_results.txt" -outfmt "6 pident qseqid sseqid evalue"
done

mkdir blast_out_serotype3a
for fasta_file in serotype_3a/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}"

    blastn -query "$fasta_file" -subject GCF_000168595.2_ASM16859v2_genomic.fna \
        -out "blast_out_serotype3a/${filename_no_ext}_serotype3a_blast_results.txt" -outfmt "6 pident qseqid sseqid evalue"
done

mkdir blast_out_serotype4b
for fasta_file in serotype_4b/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}"

    blastn -query "$fasta_file" -subject GCF_000008285.1_ASM828v1_genomic.fna \
        -out "blast_out_serotype4b/${filename_no_ext}_serotype4b_blast_results.txt" -outfmt "6 pident qseqid sseqid evalue"
done

mkdir blast_out_serotype4d
for fasta_file in serotype_4d/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}"

    blastn -query "$fasta_file" -subject GCF_000307025.1_ASM30702v1_genomic.fna \
        -out "blast_out_serotype4d/${filename_no_ext}_serotype4d_blast_results.txt" -outfmt "6 pident qseqid sseqid evalue"
done

mkdir blast_out_serotype4e
for fasta_file in serotype_4e/*.fa; do
    filename=$(basename -- "$fasta_file")
    filename_no_ext="${filename%.*}"

    blastn -query "$fasta_file" -subject GCF_000307615.1_ASM30761v1_genomic.fna \
        -out "blast_out_serotype4e/${filename_no_ext}_serotype4e_blast_results.txt" -outfmt "6 pident qseqid sseqid evalue"
done

cd ..
conda deactivate

#CHECKM COMMANDS

mkdir -pv ./checkm/{asm,db}


"$input_directory"/* ./checkm/asm/
cd ./checkm/asm

conda activate checkm
cd ../db
wget https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz
tar zxvf checkm_data_2015_01_16.tar.gz
echo 'export CHECKM_DATA_PATH=./checkm/db' >> ~/.bashrc
source ~/.bashrc
echo "${CHECKM_DATA_PATH}"

checkm taxon_set species "Listeria monocytogenes" Lm.markers

checkm \
  analyze \
  -t 128 \
  Lm.markers \
  ./asm \
  analyze_output
checkm \
  qa \
  -f checkm.tax.qa.out \
  -o 1 \
  Lm.markers \
  -t 128 \
  analyze_output
sed 's/ \+ /\t/g' checkm.tax.qa.out > checkm.tax.qa.out.tsv
cut -f 2- checkm.tax.qa.out.tsv > tmp.tab && mv tmp.tab checkm.tax.qa.out.tsv
sed -i '1d; 3d; $d' checkm.tax.qa.out.tsv
column -ts $'\t' checkm.tax.qa.out.tsv | less -S
cd ../..
conda deactivate

#QUAST COMMANDS
conda activate qual_eval

# Set the output directory
output_dir="./quast_results"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

input_dir=$input_directory
input_dir=filtered_skesa_asm

# Run Quast for each assembly in the tool's directory
for i in $(ls -A "$input_dir"); do
	mkdir -p "$output_dir/$i"
	quast.py -o "$output_dir/$i" "$input_dir/$i"/filtered_contigs.fa -r GCF_000196035.1_ASM19603v1_genomic.fna 1> "$output_dir/$i"/quast.stdout.txt 2> "$output_dir/$i"/quast.stdout.err
done

# Set the output directory
output_dir="./quast_results_3"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

input_dir=filtered_skesa_asm_3

# Run Quast for each assembly in the tool's directory
for i in $(ls -A "$input_dir"); do
	mkdir -p "$output_dir/$i"
	quast.py -o "$output_dir/$i" "$input_dir/$i"/filtered_contigs.fa -r GCF_000307025.1_ASM30702v1_genomic.fna 1> "$output_dir/$i"/quast.stdout.txt 2> "$output_dir/$i"/quast.stdout.err
done


# Deactivate the 'qual_eval' environment
conda deactivate

