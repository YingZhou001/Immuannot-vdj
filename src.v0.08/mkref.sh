
DAT=$1
PREF=$2

scriptpath=$(dirname $0)

#download_link=https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z
#dat=dat.Z
#wget ${download_link} -O ${dat}

date_tag=$(stat ${DAT} | grep Modify | cut -f2 -d' ')
refdir=${PREF}.${date_tag}
blacklist=${scriptpath}/black.list

mkdir -p ${refdir}
# prepare ref for assembly annotation



zcat ${DAT} | python3 ${scriptpath}/fmt-vdjc-gene.py "complete,c" ${blacklist}\
  | python3 ${scriptpath}/fmt-rm-dup.py \
  | gzip -c > ${refdir}/all.c.complete.fa.gz

zcat ${DAT} | python3 ${scriptpath}/fmt-vdjc-gene.py "human,complete,d,j" ${blacklist}\
  | python3 ${scriptpath}/fmt-rm-dup.py \
  | gzip -c > ${refdir}/human.dj.complete.fa.gz

zcat ${DAT} | python3 ${scriptpath}/fmt-vdjc-gene.py "human,v" ${blacklist}\
  | python3 ${scriptpath}/fmt-rm-dup.py \
  | gzip -c > ${refdir}/human.v.fa.gz

zcat ${refdir}/human.dj.complete.fa.gz ${refdir}/human.v.fa.gz \
  | python3 ${scriptpath}/fmt-sel-allele.py \
  | gzip -c > ${refdir}/human-combined-vdj.fa.gz

zcat ${refdir}/all.c.complete.fa.gz ${refdir}/human.dj.complete.fa.gz \
  ${refdir}/human.v.fa.gz \
  | python3 ${scriptpath}/fmt-sel-allele.py \
  | python3 ${scriptpath}/fmt-rm-intron.py \
  | python3 ${scriptpath}/fmt-rm-dup.py \
  | gzip -c > ${refdir}/human-combined-vdjc.nointron.fa.gz

