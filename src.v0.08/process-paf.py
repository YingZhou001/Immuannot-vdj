import sys, gzip, re
import pprint as pp


def overlap(a,b, c,d) :
    # cmpare [a,b], and [c,d]
    # assuming a<b, c<d
    if b < c or a > d :
        return(0)
    else :
        x = min(b,d) - max(a, c)
        l = min(b-a, d-c)
        return(float(x)/l)

def is_constant_gene(gene) :
    c_genes = ["IGHA","IGHA1","IGHA2","IGHA4","IGHA5","IGHD","IGHD1","IGHDD","IGHDD1P","IGHDD2P","IGHDD3P","IGHE","IGHEP1","IGHG","IGHG1","IGHG1A","IGHG1B","IGHG2","IGHG2A","IGHG2B","IGHG2C","IGHG3","IGHG3A","IGHG3B","IGHG3C","IGHG4","IGHG4A","IGHG5","IGHG5-1","IGHG5-2","IGHG6","IGHG6-1","IGHG6-2","IGHG7","IGHGP","IGHM","IGHM1","IGHM2","IGHM6","IGHMD","IGHO","IGHT1","IGHT1D","IGHT2","IGHT2D","IGHT4","IGHT5","IGHZ","IGIC1","IGIC2","IGIC3","IGIC4","IGIC5","IGKC","IGKC1","IGKC2","IGLC","IGLC1","IGLC10","IGLC11","IGLC2","IGLC2A","IGLC3","IGLC4","IGLC4P","IGLC5","IGLC6","IGLC7","IGLC8","IGLC9","IGLC/OR17-1","IGLC/OR17-2","IGLC/OR17-3","IGLC/OR22-1","IGLC/OR22-2","TRAC","TRBC1","TRBC2","TRBC3","TRDC","TRGC","TRGC1","TRGC2","TRGC2S1","TRGC3","TRGC4","TRGC5","TRGC6","TRGC7","TRGC8"]
    ret = gene.replace('"', '') in c_genes
    return(ret)



def parsing_paf(line, da_type = 'dna') :
    ds_regex = re.compile(r'ds:Z:(.*?)\s')
    llst = line.split()
    ret = {}
    ret['allele'] = llst[0]
    gene = ret['allele'].split('*')[0]
    if is_constant_gene(gene) : tag = 'C'
    else : 
        tag = gene[3]
    ret['allele_len'] = int(llst[1])
    ret['allele_mlen'] = int(llst[3]) - int(llst[2])
    reach_left_end = llst[2] == '0'
    reach_right_end = llst[3] == llst[1]
    ret['strand'] = llst[4]
    ret['ctg'] = llst[5]
    ret['ctg_len'] = int(llst[6])
    ret['ctg_fro'] = int(llst[7])
    ret['ctg_to'] = int(llst[8])
    ret['nm'] = int(llst[12].replace('NM:i:', ''))
    ret['ratio_nm'] = float(ret['nm'])/ret['allele_mlen']
    ret['pafline'] = line[:-1]
    ds_ret = ds_regex.search(line)
    if ds_ret :
        ret['ds'] = ds_ret.group(1)
        ret['cs'] = re.sub('[\[\]]', '',ret['ds'])
    else :
        ret['ds'] = ''
        ret['cs'] = ''
    # add filter here
    nmrate = float(ret['nm'])/ret['allele_mlen']
    f1 = (reach_left_end or reach_right_end)
    f1 = f1 and ret['allele_mlen']/ret['allele_len'] > 0.05
    f2 = tag == 'D' and nmrate < 0.1 and ret['allele_mlen'] > 15
    f3 = tag == 'J' and nmrate < 0.1 and ret['allele_mlen'] > 20
    if da_type == 'rna' :
        f4 = tag == 'V' and nmrate < 0.3
        if 'PART1' in ret['allele'] :
            f4 = f4 and ret['allele_mlen'] > 20
        else :
            f4 = f4 and ret['allele_mlen'] > 70
        f5 = tag == 'C' and nmrate < 0.05 and ret['allele_mlen'] > 20
    else :
        f4 = tag == 'V' and nmrate < 0.1 and ret['allele_mlen'] > 70
        f5 = tag == 'C' and nmrate < 0.1 and ret['allele_mlen'] > 50
    if f1 and (f2 or f3 or f4 or f5) :
        print('# ' + line[:-1])
        return(ret)
    else :
        return({})


def parsing_cs(cs) :
    matches = re.findall(r'([:\+\-\*]{1})([\da-z]+)', cs)
    out = [[m[0], m[1]] for m in matches]
    return(out)

def extract_cs(fro, to, cs) :
    # 0-based coordinate
    all_tab = parsing_cs(cs)
    new_tab = []
    cur_pos = 0
    cgr_l = 0
    for x in all_tab :
        tag = x[0]
        info = x[1]
        if tag  == ':' : inc = int(info)
        elif tag == '*' : inc = 1
        elif tag == '+' : inc = 0
        elif tag == '-' : inc = len(info)
        else :
            exit("unknow cigar tag: " + tag)
        cgr_l += inc
        new_pos = cur_pos + inc
        # check orverlaping between [fro, to) and [cur_pos, new_pos)
        if cur_pos >= to : break
        elif new_pos <= fro : pass
        else :
            new_l = min(to, new_pos) - max(fro, cur_pos)
            if tag == ':' : info = str(new_l)
            elif tag == '*' : info = x[1]
            elif tag == '+' : info = x[1]
            elif tag == '-' : info = x[1][0:new_l]
            new_tab.append(tag + info)
        #print("->", fro, to, cur_pos, new_pos, x, new_tab)
        cur_pos = new_pos
    return(''.join(new_tab))

def cal_mismatch(cs) :
    mat = parsing_cs(cs)
    score = 0
    for x in mat :
        tag = x[0]
        info = x[1]
        if tag == ':' : score += 0
        elif tag == '-' : score += len(info)
        elif tag == '+' : score += len(info)
        elif tag == '*' : score += 1
        else :
            exit("unknow cs tag: " + tag)
    return(score)

def compare_mapping(paf_rec1, paf_rec2) :
    map_rate1 = paf_rec1['allele_mlen']/paf_rec1['allele_len']
    map_rate2 = paf_rec2['allele_mlen']/paf_rec2['allele_len']
    if map_rate1 - map_rate2 > 0.1 : return( paf_rec1 )
    elif map_rate2 - map_rate1 > 0.1 : return( paf_rec2 )
    else : pass
    fro = max(paf_rec1['ctg_fro'], paf_rec2['ctg_fro'])
    to = min(paf_rec1['ctg_to'], paf_rec2['ctg_to'])
    mlen1 = paf_rec1['allele_mlen']
    mlen2 = paf_rec2['allele_mlen']
    if to <= fro :
        # if not overlapped, use the longer one
        if mlen1 > mlen2 : return( paf_rec1 )
        else : return (paf_rec2)
    else :
        fro1 = fro - paf_rec1['ctg_fro']
        to1 = to - paf_rec1['ctg_fro']
        subcs1 = extract_cs(fro1, to1, paf_rec1['cs'])
        fro2 = fro - paf_rec2['ctg_fro']
        to2 = to - paf_rec2['ctg_fro']
        subcs2 = extract_cs(fro2, to2, paf_rec2['cs'])
        mismatch1 = cal_mismatch(subcs1)
        mismatch2 = cal_mismatch(subcs2)
        if mismatch1 == mismatch2 : 
            if mlen1 > mlen2 : return( paf_rec1 )
            else : return (paf_rec2)
        elif mismatch1 < mismatch2 :
            return( paf_rec1 )
        else : return (paf_rec2)

def best_paf_rec(paf) :
    best_rec = paf[0]
    for x in paf :
        best_rec = compare_mapping(best_rec, x)
    return (best_rec)

pafdct = {} # define as pafdct[ctg]=[[], [], []]

da_tag = sys.argv[1]

for line in sys.stdin:
    if not line : break
    retpaf = parsing_paf(line, da_tag)
    if not retpaf : continue
    ctg = retpaf['ctg']
    if ctg not in pafdct :
        pafdct[ctg] = []
    if retpaf: 
        pafdct[ctg].append(retpaf)


# define clusters

clustdct = {}
for ctg in pafdct :
    if ctg not in clustdct :
        clustdct[ctg] = []
    sorted_pafdct = sorted(pafdct[ctg], key=lambda x: x['allele_mlen'], reverse=True)
    sorted_pafdct = sorted(pafdct[ctg], key=lambda x: x['ctg_fro'])
    for paf in sorted_pafdct :
        fro = paf['ctg_fro']
        to = paf['ctg_to']
        if not clustdct[ctg] :
            clustdct[ctg].append([fro, to])
        else :
            L = len(clustdct[ctg])
            overlap_tag = 0
            for i in range(L) :
                old_fro, old_to = clustdct[ctg][i]
                tmp = overlap(fro, to, old_fro, old_to)
                # overlaping region should larger than 50%
                if overlap(fro, to, old_fro, old_to) > 0.5:
                    clustdct[ctg][i] = [min(fro, old_fro), max(to, old_to)]
                    overlap_tag = 1
                    break
            if overlap_tag == 0:
                clustdct[ctg].append([fro, to])

# search gene pattern

for ctg in clustdct :
    output = []
    for clust in clustdct[ctg] :
        newpaf = []
        for paf in pafdct[ctg] :
            fro = paf['ctg_fro']
            to = paf['ctg_to']
            if overlap(fro, to, clust[0], clust[1]) > 0.9:
                newpaf.append(paf)
        # search the best allele for each cluster
        ## sorted paf
        if not newpaf : continue
        newpaf = sorted(newpaf, key=lambda x: x['allele_mlen'], reverse=True)

        tmp_L = len(newpaf)
        if tmp_L > 30 :
            tmp_L = int(len(newpaf)/2) # jsut use top 50% long ones to infer 
        selpaf = best_paf_rec(newpaf[0:tmp_L])

        nm_cut = selpaf['nm']
        mlen_cut = selpaf['allele_mlen']
        sel_alleles = []
        for x in newpaf :
            if x['nm'] <= nm_cut and x['allele_mlen'] >= mlen_cut :
                sel_alleles.append(x['strand']+x['allele'])

        output.append([';'.join(sel_alleles), selpaf['ctg_fro'],
                       selpaf['ctg_to'], str(selpaf['ctg_len']),
                       str(round(selpaf['allele_mlen']*100/selpaf['allele_len'],1)), 
                       str(selpaf['nm']), selpaf['pafline']])
    output = sorted(output, key=lambda x: int(x[2]))
    for x in output :
        ostr=ctg + '\t' + str(x[1]) + '\t' +  str(x[2])
        ostr+= '\t' + x[0] + '\t' + x[3]
        ostr+= '\tmapping_rate=' + x[4] + "%"
        ostr+= ',nm=' + x[5]
        ostr = ostr.replace('"', '')
        ostr+= '\t@' + x[6]
        print(ostr)


