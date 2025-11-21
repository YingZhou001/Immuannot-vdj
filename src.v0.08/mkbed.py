import sys

gap_cut = 1000000
extend = 50000

ctgs = {}
length = {}

for line in sys.stdin:
    if not line: break
    fields = line.strip().split("\t")
    allele_id = fields[0]
    allele_len = int(fields[1])
    ctg = fields[5]
    ctg_len = int(fields[6])
    ctg_fro = int(fields[7])
    ctg_to = int(fields[8])
    mlen = int(fields[10])
    reg = allele_id[0:3]
    length[ctg] = ctg_len
    if mlen/allele_len > 0.3 :
        if ctg not in ctgs : ctgs[ctg] = []
        ctgs[ctg].append([ctg_fro, ctg_to, [reg]])

def merge_interval(x, gap_cut) :
    arr = sorted(x, key = lambda x: x[0])
    arr_l = len(arr)
    new_arr = [arr[0]]
    for i in range(1, arr_l) :
        fro,to,reg = arr[i]
        old_fro,old_to,old_reg = new_arr[-1]
        if old_to + gap_cut > fro :
            new_arr[-1][1] = to
            if reg[0] not in old_reg :
                new_arr[-1][2].append(reg[0])
        else :
            new_arr.append(arr[i])
    return(new_arr)

for ctg in ctgs :
    arr = merge_interval(ctgs[ctg], gap_cut)
    for a,b,c in arr :
        print(ctg, max(0, a - extend), min(b + extend, length[ctg]), ','.join(c), sep='\t')
