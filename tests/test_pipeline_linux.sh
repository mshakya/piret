#!/usr/bin/env bash

blue=$(tput setaf 4)
normal=$(tput sgr0)

if [[ ( $1 -eq "prok") ]];
then
	
	# create a test experimental design file for proks
    printf "${blue}Creating experimental design file for prokaryotes\n${normal}"
    find tests/data/fastqs -iname "*test*R2.fastq" | sort > testR2.txt
    find tests/data/fastqs -iname "*test*R1.fastq" | sort | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
    printf "samp1\nsamp2\nsamp3\nsamp4\nsamp5\nsamp6\n" > test_id.txt
    printf "liver\nspleen\nspleen\nliver\nliver\nspleen\n" > test_gr.txt
    paste test_id.txt testR1.txt testR2.txt test_gr.txt > test_ed.txt
    rm test_id.txt testR1.txt testR2.txt test_gr.txt
    sed "s/:$(printf '\t')/:/g" test_ed.txt > test_edII.txt
    printf "#SampleID\tFiles\tGroup\n" | cat - test_edII.txt > test_prok.txt
    rm test_ed.txt test_edII.txt

fi

if [[ ( $1 -eq "prok" && $2 -eq "hisat" ) || ( -z $1 ) ]];
then
    #==============================================================================#
    printf "${blue}testing pipeline with prokarya using DEedge\n${normal}"

    printf "bin/piret -d tests/test_prok -e test_prok.txt \
    -gp tests/data/test_prok.gff \
    -i tests/test_prok/prok_index -k prokarya -m DEedge \
    -fp tests/data/test_prok.fna --config luigi.cfg"

    bin/piret -d tests/test_prok -e test_prok.txt \
    -gp tests/data/test_prok.gff \
    -i tests/test_prok/prok_index -k prokarya -m DEedge \
    -fp tests/data/test_prok.fna --config luigi.cfg \
    --aligner hisat2


    printf "${blue}fininshed testing pipeline for prokarya only\n${normal}"
    #==============================================================================#
fi
	############################################################################


# euk star ballgown
bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 -k eukarya -m ballgown -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa --aligner STAR --config luigi.cfg 


# euk star novel ballgown
bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 -k eukarya -m ballgown -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa --aligner STAR --config luigi.cfg --novel

# euk star edgeR
 bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 -k eukarya -m edgeR -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa --aligner STAR --config luigi.cfg 

# euk star edger novel (fails at pathway step)

# euk star deseq

# euk hisat2 edgeR novel
bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 -k eukarya -m edgeR -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa --aligner hisat2 --config luigi.cfg --novel

# euk hisat2 edgeR
bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 -k eukarya -m edgeR -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa --aligner hisat2 --config luigi.cfg

# both hisat

if [[ $1 -eq "euk" ]] || [[ -z $1 ]];
then
    printf "${blue}testing pipeline with eukarya \n${normal}"

    find tests/data2/fastqs -iname "*read2.fastq" | sort > testR2.txt
    find tests/data2/fastqs -iname "*read1.fastq" | sort | sed 's/read1\.fastq/read1\.fastq:/g' > testR1.txt
    printf "samp1\nsamp2\nsamp3\nsamp4\nsamp5\nsamp6\n" > test_id.txt
    printf "HBR\nHBR\nHBR\nUBR\nUBR\nUBR\n" > test_gr.txt
    paste test_id.txt testR1.txt testR2.txt test_gr.txt > test_ed.txt
    rm test_id.txt testR1.txt testR2.txt test_gr.txt
    sed "s/:$(printf '\t')/:/g" test_ed.txt > test_edII.txt
    printf "#SampleID\tFiles\tGroup\n" | cat - test_edII.txt > test_euk.txt
    rm test_ed.txt test_edII.txt
fi


# printf "bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 \
# -k eukarya -m DESeq2 -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa \
# --novel --aligner STAR --config luigi.cfg "

# bin/piret -d tests/test_euk -e test_euk.txt -ge tests/data2/chr22_ERCC92.gff3 \
# -k eukarya -m DESeq2 -fe tests/data2/chr22_ERCC92.fa  -p 0.05 --org hsa \
# --novel --aligner STAR --config luigi.cfg 

# bin/piret -d tests/test_euk -e test_euk.txt \
# -ge tests/data2/chr22_ERCC92.gff3 \
# -i tests/test_euk/euk_index -k eukarya -m all \
# -fe tests/data2/chr22_ERCC92.fa --config luigi.cfg -p 0.05

# printf "${blue}fininshed testing pipeline for eukarya\n${normal}"

# # #==============================================================================#
# # printf "${blue}testing pipeline with prok and euk now \n${normal}"

# printf "bin/piret -d tests/test_both -e test_prok.txt \
# -gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
# -i tests/test_both/prok_euk_index -k both -m all \
# -fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa"

# # mkdir -p tests/test_both/samp1/trimming_results 
# # mkdir -p tests/test_both/samp2/trimming_results 
# # mkdir -p tests/test_both/samp3/trimming_results 
# # mkdir -p tests/test_both/samp4/trimming_results 
# # mkdir -p tests/test_both/samp5/trimming_results 
# # mkdir -p tests/test_both/samp6/trimming_results 

# # cp tests/test_both/samp1/trimming_results/*.* tests/test_euk/samp1/trimming_results/
# # cp tests/test_both/samp2/trimming_results/*.* tests/test_euk/samp2/trimming_results/ 
# # cp tests/test_both/samp3/trimming_results/*.* tests/test_euk/samp3/trimming_results/ 
# # cp tests/test_both/samp4/trimming_results/*.* tests/test_euk/samp4/trimming_results/ 
# # cp tests/test_both/samp5/trimming_results/*.* tests/test_euk/samp5/trimming_results/ 
# # cp tests/test_both/samp6/trimming_results/*.* tests/test_euk/samp6/trimming_results/ 

bin/piret -d tests/test_both -e test_prok.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m all \
-fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa

# printf "${blue}fininshed testing pipeline for prok and euk\n${normal}"


# # printf "${blue}Cleanning up!!\n${normal}"
# # # rm -rf test_experimental_design.txt
# # # rm -rf tests/test_both tests/test_euk tests/test_prok

printf "${green}piret with both euk and prok, and edgeR\n${normal}"
bin/piret -d tests/test_both_hisat2_edger -e test_prok.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m edgeR \
-fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa --novel
printf "${blue}fininshed testing pipeline for prok and euk\n${normal}"


bin/piret -d tests/test_both_hisat2_edger -e test_prok.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m edgeR \
-fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa --novel \
--aligner STAR

bin/piret -d tests/test_both_hisat2_ballgown -e test_prok.txt \
-gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
-i tests/test_both/prok_euk_index -k both -m ballgown \
-fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa --novel

 bin/piret -d tests/test_both_hisat2_deseq2 -e test_prok.txt \
 -gp tests/data/test_prok.gff -ge tests/data/eukarya_test.gff3 \
 -i tests/test_both/prok_euk_index -k both -m DESeq2 \
 -fp tests/data/test_prok.fna -fe tests/data/eukarya_test.fa --novel