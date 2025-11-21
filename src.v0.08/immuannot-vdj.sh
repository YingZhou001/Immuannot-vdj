#!/bin/bash - 
#===============================================================================
#
#          FILE: immuannot-vdj.sh
# 
#         USAGE: bash ./immuannot-vdj.sh  [OPTIONS] value
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Ying Zhou
#  ORGANIZATION: 
#       CREATED: 05/20/2025 15:32
#      REVISION: v0.08
#===============================================================================

set -o nounset                              # Treat unset variables as an error

SCRIPTPATH=$(dirname $0)

MINIMAP2=${SCRIPTPATH}/dep-tools/bin/minimap2
SEQTK=${SCRIPTPATH}/dep-tools/bin/seqtk

MODE=unset
REFDIR=unset
TARGET=unset
OUTDIR=./output
THREAD=3




usage()
{
  echo "
  Usage: bash ${SCRIPTPATH}/immuannot-vdj.sh  [OPTIONS] value
      [ --asm/--hifidna/--rna/--ontuldna   flag indicating assembly or hifi reads as input     ]
      [ -r | --refdir <string> references directory                             ]
      [ -t | --target <string> target input (.fa, fa.gz, .fq, .fq.gz)           ]
      [ -o | --outdir <string> output directory (optional, default ./output)    ]
      [ -n | --nthread <int>   num of thread (optional, default 3)              ]
  "
  exit 2
}


PARSED_ARGUMENTS=$(getopt -a -n immuannot-vdj \
  -o r:t:o:n: \
  --long refdir:,target:,outdir:,nthread:,asm,hifidna,rna,ontuldna\
  -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
  exit 1
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -r | --refdir)  REFDIR="$2"    ; shift 2 ;;
    -t | --target)  TARGET="$2"    ; shift 2 ;;
    -o | --outdir)  OUTDIR="$2"    ; shift 2 ;;
    -n | --nthread) THREAD="$2"    ; shift 2 ;;
    --asm)  MODE="ASM"    ; shift ;;
    --hifidna)  MODE="hifiDNA"    ; shift ;;
    --ontuldna)  MODE="ONT-UL"    ; shift ;;
    --rna)  MODE="RNA"    ; shift ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1."
      usage
      ;;
  esac
done


mkdir -p ${OUTDIR}

if [ $MODE == unset ]; then
  echo "Error: '--asm/--hifidna/--ontuldna/--rna' is required."
  usage
fi

if [ $REFDIR == unset ]; then
  echo "Error: '-r/--refdir' is required."
  usage
fi

if [ $TARGET == unset ]; then
  echo "Error: '-t/--target' is required."
  usage
fi


start_second=`date +%s`
start=`date +%D-%H:%M:%S`

echo "########################################"
echo "##Welcome###############################"
echo "########################################"
echo ""
echo "Starting time: ${start}"
echo "#####parameters:########################"
echo "TARGET(-t)                       : $TARGET"
echo "MODE(--asm/--hifidna/--ontuldna/--rna)      : $MODE"
echo "REFDIR(-r)                       : $REFDIR"
echo "OUTDIR(-o)                       : $OUTDIR"
echo "nTHREAD(-n)                      : $THREAD"


pref=${OUTDIR}/tmp
out_paf=${pref}.mm2.paf.gz
annot0=${pref}.annot.gz

ref_dna_vdj=${REFDIR}/human-combined-vdj.fa.gz
ref_dna_c=${REFDIR}/all.c.complete.fa.gz
ref_rna_vdjc=${REFDIR}/human-combined-vdjc.nointron.fa.gz

#dev-options
run_extraction=true
run_aligment=true
run_annot=true



if ${run_extraction}
then
  if [ ${MODE} == 'ASM' ]
  then
    echo "extract regions from assembly to run faster"
    bed=${pref}.bed
    TARGET_fa=${pref}.fa.gz
    (${MINIMAP2} -t${THREAD} ${TARGET} ${ref_dna_vdj}; \
      ${MINIMAP2} -t${THREAD} ${TARGET} ${ref_dna_c}) \
      | python3 ${SCRIPTPATH}/mkbed.py > ${bed}
    ${SEQTK} subseq ${TARGET} ${bed} | gzip -c > ${TARGET_fa}
    TARGET=${TARGET_fa}
  elif [ ${MODE} == 'ONT-UL' ]
  then
    echo "extract regions from assembly to run faster"
    bed=${pref}.bed
    TARGET_fa=${pref}.fa.gz
    (${MINIMAP2} -t${THREAD} ${TARGET} ${ref_dna_vdj}; \
      ${MINIMAP2} -t${THREAD} ${TARGET} ${ref_dna_c}) \
      | python3 ${SCRIPTPATH}/mkbed.py > ${bed}
    ${SEQTK} subseq ${TARGET} ${bed} | ${SEQTK} seq -A | gzip -c > ${TARGET_fa}
    TARGET=${TARGET_fa}
  fi
fi

if ${run_aligment}
then
  options1="-t${THREAD} -c -k9 -I10 -n2 -w3 -s25 -m13 --end-bonus=5 -f1000 -N100 -z50,25 --ds -v2"
  options2="-t${THREAD} -c -k19 -I10 -n2 -w19 -s100 -m13 --end-bonus=5 -f1000 -N100 -z50,25 --ds -v2"
  options3="-t${THREAD} -cx map-ont -I10 --ds"
  if [ ${MODE} == 'ASM' ]
  then
    (${MINIMAP2} ${options1} ${TARGET} ${ref_dna_vdj};  \
      ${MINIMAP2} ${options2} ${TARGET} ${ref_dna_c}) \
      | python3 ${SCRIPTPATH}/process-paf.py 'dna' | gzip -c > ${annot0}
  elif [ ${MODE} == 'hifiDNA' ]
  then
    (${MINIMAP2} ${options1} ${TARGET} ${ref_dna_vdj};  \
      ${MINIMAP2} ${options2} ${TARGET} ${ref_dna_c}) \
      | python3 ${SCRIPTPATH}/process-paf.py 'dna' | gzip -c > ${annot0}
  elif [ ${MODE} == 'ONT-UL' ]
  then
    (${MINIMAP2} ${options1} ${TARGET} ${ref_dna_vdj};  \
      ${MINIMAP2} ${options3} ${TARGET} ${ref_dna_c}) \
      | python3 ${SCRIPTPATH}/process-paf.py 'dna' | gzip -c > ${annot0}
 
  elif [ ${MODE} == 'RNA' ]
  then
    ${MINIMAP2} ${options1} ${TARGET} ${ref_rna_vdjc} \
      | python3 ${SCRIPTPATH}/process-paf.py 'rna' | gzip -c > ${annot0}
  fi
fi


#annotation
if ${run_annot}
then
  pref=${OUTDIR}/${MODE}.annot
  out_gff=${pref}.gff.gz
  hap_summary=${pref}.summary.txt
  info=${pref}.info.gz

  if [[ ${MODE} == 'RNA' ]]
  then
    zgrep '>' ${ref_rna_vdjc} | gzip -c > ${info}
    zgrep -v '^#' ${annot0} \
      | python3 ${SCRIPTPATH}/add-elements.py ${info} \
      | gzip -c > ${out_gff}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/read-smry-rna.py > ${hap_summary}
  elif [[ ${MODE} == 'hifiDNA' ]] 
  then
    zgrep '>' ${ref_dna_vdj} | gzip -c > ${info}
    zgrep '>' ${ref_dna_c} | gzip -c >> ${info}
    zgrep -v '^#' ${annot0} \
      | python3 ${SCRIPTPATH}/add-elements.py ${info} \
      | gzip -c > ${out_gff}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/read-smry-dna.py > ${hap_summary}
  elif [[ ${MODE} == 'ONT-UL' ]] 
  then
    zgrep '>' ${ref_dna_vdj} | gzip -c > ${info}
    zgrep '>' ${ref_dna_c} | gzip -c >> ${info}
    zgrep -v '^#' ${annot0} \
      | python3 ${SCRIPTPATH}/add-elements.py ${info} \
      | gzip -c > ${out_gff}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/read-smry-dna.py > ${hap_summary}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/asm-gff-restore.py | gzip -c > ${pref}.tmp
    mv ${pref}.tmp ${out_gff}
  elif [[ ${MODE} == 'ASM' ]] 
  then
    zgrep '>' ${ref_dna_vdj} | gzip -c > ${info}
    zgrep '>' ${ref_dna_c} | gzip -c >> ${info}
    zgrep -v '^#' ${annot0} \
      | python3 ${SCRIPTPATH}/add-elements.py ${info} \
      | gzip -c > ${out_gff}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/read-smry-dna.py > ${hap_summary}
    zcat ${out_gff} | python3 ${SCRIPTPATH}/asm-gff-restore.py | gzip -c > ${pref}.tmp
    mv ${pref}.tmp ${out_gff}
  fi
  rm ${info}
fi

end_second=`date +%s`
end=`date +%D-%H:%M:%S`
runtime=$((end_second - start_second))
echo "Ending time: ${end}, Wallclock time :${runtime} seconds"
echo "####done####"
