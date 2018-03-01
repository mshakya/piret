#!/usr/bin/env bash

blue=$(tput setaf 4)
normal=$(tput sgr0)
# create a test experimental design file

# find tests/data/fastqs -iname "*test*R2.fastq" -exec readlink -f {} \; | sort > testR2.txt
# find tests/data/fastqs -iname "*test*R1.fastq" -exec readlink -f {} \; | sort | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
find tests/data/fastqs -iname "*test*R2.fastq" | sort > testR2.txt
find tests/data/fastqs -iname "*test*R1.fastq" | sort | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
printf "samp1\nsamp2\nsamp3\nsamp4\nsamp5\nsamp6\n" > test_id.txt
printf "liver\nspleen\nspleen\nliver\nliver\nspleen\n" > test_gr.txt
paste test_id.txt testR1.txt testR2.txt test_gr.txt > test_ed.txt
rm test_id.txt testR1.txt testR2.txt test_gr.txt
sed -i 's/:\t/:/g' test_ed.txt
printf "ID\tFiles\tGroup\n" | cat - test_ed.txt > test_experimental_design.txt
rm test_ed.txt

#==============================================================================#
printf "${blue}testing pipeline with prokarya only\n${normal}"

printf "bin/runPiReT -d tests/test_prok -e test_experimental_design.txt \
-gp tests/data/test_prok.gff \
-i tests/test_prok/prok_index -k prokarya -m all \
-fp tests/data/test_prok.fa"

bin/runPiReT -d tests/test_prok -e test_experimental_design.txt \
-gp tests/data/test_prok.gff \
-i tests/test_prok/prok_index -k prokarya -m all \
-fp tests/data/test_prok.fa

printf "${blue}fininshed testing pipeline for prokarya only\n${normal}"
#==============================================================================#
printf "${blue}testing pipeline with eukarya now\n${normal}"

printf "bin/runPiReT -d tests/test_euk -e test_experimental_design.txt \
-ge tests/data/eukarya_test.gff3 \
-i tests/test_euk/euk_index -k eukarya -m all \
-fe tests/data/eukarya_test.fa"

mkdir -p tests/test_euk/samp1/trimming_results 
mkdir -p tests/test_euk/samp2/trimming_results 
mkdir -p tests/test_euk/samp3/trimming_results 
mkdir -p tests/test_euk/samp4/trimming_results 
mkdir -p tests/test_euk/samp5/trimming_results 
mkdir -p tests/test_euk/samp6/trimming_results 

cp tests/test_prok/samp1/trimming_results/*.* tests/test_euk/samp1/trimming_results/
cp tests/test_prok/samp2/trimming_results/*.* tests/test_euk/samp2/trimming_results/ 
cp tests/test_prok/samp3/trimming_results/*.* tests/test_euk/samp3/trimming_results/ 
cp tests/test_prok/samp4/trimming_results/*.* tests/test_euk/samp4/trimming_results/ 
cp tests/test_prok/samp5/trimming_results/*.* tests/test_euk/samp5/trimming_results/ 
cp tests/test_prok/samp6/trimming_results/*.* tests/test_euk/samp6/trimming_results/ 

bin/runPiReT -d tests/test_euk -e test_experimental_design.txt \
-ge tests/data/eukarya_test.gff3 \
-i tests/test_euk/euk_index -k eukarya -m all \
-fe tests/data/eukarya_test.fa

printf "${blue}fininshed testing pipeline for eukarya\n${normal}"

#==============================================================================#
printf "${blue}testing pipeline with prok and euk now \n${normal}"

printf "bin/runPiReT -d tests/test_both -e test_experimental_design.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m all \
-fp tests/data/test_prok.fa -fe tests/data/eukarya_test.fa\n"

mkdir -p tests/test_both/samp1/trimming_results 
mkdir -p tests/test_both/samp2/trimming_results 
mkdir -p tests/test_both/samp3/trimming_results 
mkdir -p tests/test_both/samp4/trimming_results 
mkdir -p tests/test_both/samp5/trimming_results 
mkdir -p tests/test_both/samp6/trimming_results 

cp tests/test_both/samp1/trimming_results/*.* tests/test_euk/samp1/trimming_results/
cp tests/test_both/samp2/trimming_results/*.* tests/test_euk/samp2/trimming_results/ 
cp tests/test_both/samp3/trimming_results/*.* tests/test_euk/samp3/trimming_results/ 
cp tests/test_both/samp4/trimming_results/*.* tests/test_euk/samp4/trimming_results/ 
cp tests/test_both/samp5/trimming_results/*.* tests/test_euk/samp5/trimming_results/ 
cp tests/test_both/samp6/trimming_results/*.* tests/test_euk/samp6/trimming_results/ 

bin/runPiReT -d tests/test_both -e test_experimental_design.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m all \
-fp tests/data/test_prok.fa -fe tests/data/eukarya_test.fa

printf "${blue}fininshed testing pipeline for prok and euk\n${normal}"


printf "${blue}Cleanning up!!\n${normal}"
# rm -rf test_experimental_design.txt
# rm -rf tests/test_both tests/test_euk tests/test_prok