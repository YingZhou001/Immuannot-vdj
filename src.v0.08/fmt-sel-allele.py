import sys,textwrap

# usage :
# zcat inp.fa.gz | python3 sel.py > out.txt

## rules
## 1) keep one allele if alleles share the same seq and allele type
## 2) keep the longer allele if the short one is perfectly contained by the
##      longer one
## 3) keep alleles with the same allele name but different seqs

# load fa

seqs = {}
seq = ''

for line in sys.stdin:
    if not line : break
    if line[0] == '>' :
        if seq :
            if seq not in seqs: seqs[seq] = []
            seqs[seq].append(rec)
            seq = ''
        rec = line[1:-1]
    else :
        seq += line[:-1]


print("input total num of diff seqs:", len(seqs), file=sys.stderr)

seq_arr = []

for seq in seqs :
    seq_arr.append(seq)

seq_arr = sorted(seq_arr, key = lambda x : len(x), reverse = True)

L = len(seq_arr)

rmlist_short = []

for i in range(L-1) :
    seq1 = seq_arr[i]
    for j in range(i+1, L) :
        seq2 = seq_arr[j]
        if seq2 in seq1 :
            if seq2 in seqs :
                rmlist_short.append(seqs[seq2])
                seqs.pop(seq2)

rmlist_dup = []

sel_list = []

for seq in seqs :
    alleles = {}
    for x in seqs[seq] :
        allele = x.split('|')[0]
        if allele not in alleles : alleles[allele] = x
        else :
            rmlist_dup.append(x)
    for allele in alleles:
        sel_list.append(alleles[allele])


for seq in seqs :
    for x in seqs[seq] :
        if x in sel_list :
            print('>' + x)
            print(textwrap.fill(seq, width=80))

