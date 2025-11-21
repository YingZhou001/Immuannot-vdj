import sys

da = {}

key_bag=["L-PART1", "L-PART2", "V-REGION", "V-HEPTAMER", "V-SPACER", "V-NONAMER",
         "CH1","CH2","CH3","CH4","CHS", "CH5", "CH6", "CH7",
         "CH1D","CH2D","CH2D2","CH2D3","CH3D","CH3D2","CH3D3",
         "CH4D","CH4D2","CH4D3", "CH8", "CH9", "CH10","CHT","CHX",
         "EX0", "EX4", "EX5", "EX6", "EX7", "EX8", "EX9",
         "EX1","EX2","EX2A","EX2B","EX2R","EX2T","EX3","EX4UTR",
         "H","I-EXON","M","M1","M2", "H1", "H2", "H3", "H4", "H5",
         "STOP-CODON", "5'UTR", "3'UTR","C-REGION",
        ]

def load_meta(line, key_bag=key_bag) :
    ret = {}
    llst = line.split()
    ret['ctg'] = llst[0][1:]
    ret['imgt_type'] = llst[1]
    ret['ele'] = []
    for x in llst[2:] :
        k,v = x.split('=')
        if k == 'tot_length' : 
            ret['tot_length'] = v
        else :
            if k not in key_bag : continue
            rg0,rg1 = v.split('..')
            rg0,rg1 = int(rg0), int(rg1)
            if k == "3'UTR" :
                rg1 = min(rg1, rg0+50)
            elif k == "5'UTR" :
                rg0 = max(rg0, rg1-50)
            ret['ele'].append([k, int(rg0), int(rg1)])
    ret['ele'] = sorted(ret['ele'], key = lambda x: x[1])
    return(ret) 

def merge_ele(ele_arr) :
    new_arr = []
    for x,fro,to in ele_arr :
        if not new_arr : new_arr.append([x, fro, to])
        else :
            if fro <= new_arr[-1][2] + 1 :
                new_arr[-1][0] += '~' + x
                new_arr[-1][2] = max(new_arr[-1][2], to)
            else :
                new_arr.append([x, fro, to])
    ele_clust = []
    for x,fro,to in new_arr :
        ele_clust.append(x.split('~'))
    return(ele_clust)

def extract_ele(seq, ele_arr) :
    # ele_arr should be sorted
    ## extract seq
    newseq = ''
    maparr = []
    fro_old = -1
    to_old = -1
    cur_pos = 1
    for x,fro,to in ele_arr :
        fro = fro
        to = to
        if fro <= to_old : fro = to_old + 1
        if to < fro : continue
        newseq += seq[fro-1:to]
        l = to - fro + 1
        maparr.append([fro, to, cur_pos, cur_pos+l-1])
        #newarr.append([x, str(cur_pos), str(cur_pos+l-1)])
        cur_pos += l 
        fro_old = fro
        to_old = to
    ## re-coordinate arr
    newarr = []
    for x,fro,to in ele_arr :
        for a,b,c,d in maparr :
            if a <= fro and fro <= b : new_fro = c + fro - a
            if a <= to and to <= b : new_to = c + to - a
        newarr.append([x, new_fro, new_to])
    return([newseq, newarr])


def split_seq(seq, ele_arr, clust) :
    sub_ele_arr = []
    for x,fro,to in ele_arr:
        if x in clust :
            sub_ele_arr.append([x,fro,to])

    newseq, new_ele = extract_ele(seq, sub_ele_arr)
    metaline = []
    eletag = '-'.join(clust)
    alelle,source = meta['ctg'].split('|')
    newctg = alelle + '=' + eletag + '|' + source
    metaline.append('>'+newctg)
    metaline.append(meta['imgt_type'])
    for k,rg0,rg1 in new_ele :
        metaline.append(k + '=' + str(rg0) + '..' + str(rg1))
    metaline.append('tot_length' + '=' + str(len(newseq)))
    return([newseq, ' '.join(metaline)])

seq = ''


for line in sys.stdin:
    if not line: break
    if line[0] == '>' :
        if seq :
            # for D and J segment no element is extracted
            if not meta['ele'] :
                print(metaline[:-1])
                print(seq)
            else :
                ele_clust = merge_ele(meta['ele'])

                for clust in ele_clust :
                    newseq, metaline = split_seq(seq, meta['ele'], clust)
                    print(metaline)
                    print(newseq)

            seq = ''
        metaline = line
        meta = load_meta(line)
    else :
        seq += line[:-1]

# last item
if not meta['ele'] :
    print(metaline[:-1])
    print(seq)
else :
    ele_clust = merge_ele(meta['ele'])
    for clust in ele_clust :
        newseq,metaline = split_seq(seq, meta['ele'], clust)
        print(metaline)
        print(newseq)
