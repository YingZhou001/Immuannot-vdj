#!/bin/bash - 
#===============================================================================
#
#          FILE: test-dna.sh
# 
#         USAGE: ./test-dna.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 05/20/2025 13:57
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

script=../src/immuannot-vdj.sh
refdir=ref.2025-11-21

REF=true
DNA=false
RNA=false
ASM=false

test_mkref(){
	dat=dat.Z
	# wget https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z -O ${dat}
	bash ../src/mkref.sh ${dat} ref

}


test_hifiDNA(){
	## test dna reads
	tar=dna.fa.gz
	outdir=hifi_dna-test
	nthread=10
	time bash ${script} --hifidna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}
}

test_ulDNA(){
	## test dna reads
	tar=ul_dna.fa.gz
	outdir=ul_dna-test
	nthread=10
	time bash ${script} --ontuldna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}
}


test_RNA(){
	## test rna transcript
	tar=rna.fa.gz
	outdir=rna-test
	nthread=10
	time bash ${script} --rna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}
}


test_ASM(){
	## test dna reads
	tar=asm.fa.gz
	outdir=asm-test
	nthread=10
	time bash ${script} --asm -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}
}

#test_mkref

test_ulDNA
#test_hifiDNA
#test_RNA
#test_ASM
