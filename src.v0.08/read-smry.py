import sys,re

c_genes = ["IGHA","IGHA1","IGHA2","IGHA4","IGHA5","IGHD","IGHD1","IGHDD","IGHDD1P","IGHDD2P","IGHDD3P","IGHE","IGHEP1","IGHG","IGHG1","IGHG1A","IGHG1B","IGHG2","IGHG2A","IGHG2B","IGHG2C","IGHG3","IGHG3A","IGHG3B","IGHG3C","IGHG4","IGHG4A","IGHG5","IGHG5-1","IGHG5-2","IGHG6","IGHG6-1","IGHG6-2","IGHG7","IGHGP","IGHM","IGHM1","IGHM2","IGHM6","IGHMD","IGHO","IGHT1","IGHT1D","IGHT2","IGHT2D","IGHT4","IGHT5","IGHZ","IGIC1","IGIC2","IGIC3","IGIC4","IGIC5","IGKC","IGKC1","IGKC2","IGLC","IGLC1","IGLC10","IGLC11","IGLC2","IGLC2A","IGLC3","IGLC4","IGLC4P","IGLC5","IGLC6","IGLC7","IGLC8","IGLC9","IGLC/OR17-1","IGLC/OR17-2","IGLC/OR17-3","IGLC/OR22-1","IGLC/OR22-2","TRAC","TRBC1","TRBC2","TRBC3","TRDC","TRGC","TRGC1","TRGC2","TRGC2S1","TRGC3","TRGC4","TRGC5","TRGC6","TRGC7","TRGC8"]

db = {}
read_info = {}

for line in sys.stdin:
    if not line : break
    if line[0] == '#': continue
    llst = line.split('\t')
    if llst[2] == 'transcript' : continue
    attr = {}
    for x in llst[-1][:-1].split('; ') :
        k,v = x.split(' ')
        attr[k] = v
    ctg = llst[0]
    if ctg not in db : db[ctg] = {}
    gene_id = attr['gene_id']
    gene_name = attr['gene_name']
    if gene_id not in db[ctg] :
        db[ctg][gene_id] = {}
        db[ctg][gene_id][gene_name] = {}
    feature = llst[2]
    fro = llst[3]
    to = llst[4]
    strand = llst[6]
    mapping = attr['mapping_rate']
    if 'gene' in feature :
        ds_nm = attr['ds_nm']
        db[ctg][gene_id][gene_name]['gene'] = [fro, to, strand, mapping, ds_nm]
    else :
        db[ctg][gene_id][gene_name][feature] = [fro, to, strand, mapping]


def check_part(da, tag) :
    strand = da['gene'][2]
    ret = ''
    if tag == 'J' :
        is_hpt = "J-HEPTAMER" in da
        is_spc = "J-SPACER" in da
        is_nnm = "J-NONAMER" in da
        is_rss = is_hpt and is_spc and is_nnm
        is_cds = False
        if "J-REGION" in da :
            fro,to = da['J-REGION'][0:2]
            if abs(int(to) - int(fro)) > 10 : is_cds = True
        if strand == '-' :
            if is_cds : ret += 'j'
            if is_rss : ret += '>'
        else :
            if is_rss : ret += '<'
            if is_cds : ret += 'J'
    elif tag == 'D' :
        is_hpt5 = "5'D-HEPTAMER" in da
        is_spc5 = "5'D-SPACER" in da
        is_nnm5 = "5'D-NONAMER" in da
        is_5rss = is_hpt5 and is_spc5 and is_nnm5

        is_hpt3 = "3'D-HEPTAMER" in da
        is_spc3 = "3'D-SPACER" in da
        is_nnm3 = "3'D-NONAMER" in da
        is_3rss = is_hpt3 and is_spc3 and is_nnm3

        is_cds = False
        if "D-REGION" in da :
            fro,to = da['D-REGION'][0:2]
            if abs(int(to) - int(fro)) > 5 : is_cds = True
        if strand == '-' :
            if is_3rss : ret += '<'
            if is_cds : ret += 'd'
            if is_5rss : ret += '>'
        else :
            if is_5rss : ret += '<'
            if is_cds : ret += 'D'
            if is_3rss : ret += '>'

    elif tag == 'V' :
        is_hpt = "V-HEPTAMER" in da
        is_spc = "V-SPACER" in da
        is_nnm = "V-NONAMER" in da
        is_rss = is_hpt and is_spc and is_nnm

        is_ldr = "L-PART2" in da

        is_cds = False
        if "V-REGION" in da :
            fro,to = da['V-REGION'][0:2]
            if abs(int(to) - int(fro)) > 50 : is_cds = True
        
        if strand == '-' :
            if is_rss : ret += '<'
            if is_cds : ret += 'v'
            if is_ldr : ret += '$'
        else :
            if is_ldr : ret += '$'
            if is_cds : ret += 'V'
            if is_rss : ret += '>'
    return(ret)


def check_strand(da) :
    pos_l = 0
    neg_l = 0
    for gene_id in da :
        for gene_name in da[gene_id] :
            fro,to,strand = da[gene_id][gene_name]['gene'][0:3]
            if strand == '+' : pos_l += int(to) - int(fro)
            else : neg_l += int(to) - int(fro)
    if pos_l >= neg_l : return('+')
    else : return('-')

def rev_hap_str(hap_str) :
    rev_map = {'>': '<', '<':'>', 
               'j': 'J', 'J':'j',
               'v': 'V', 'V':'v',
               'd': 'D', 'D':'d',
               '$': '$'}
    new_str = []
    for x in hap_str.split('~')[::-1] :
        new_x = []
        for y in x.split('|')[::-1] :
            if y[0] not in rev_map :
                new_x.append(y)
            else :
                new_x.append(''.join([rev_map[z] for z in y[::-1]]))
        new_str.append('|'.join(new_x))
    new_str = '~'.join(new_str)
    return(new_str)

for ctg in db:
    vdj_tag = 'F'
    hap = []
    ctg_strand = check_strand(db[ctg])
    for gene_id in db[ctg] :
        for gene_name in db[ctg][gene_id] :
            tag = gene_name[4]
            ret = check_part(db[ctg][gene_id][gene_name], tag)
            fro,to = db[ctg][gene_id][gene_name]['gene'][0:2]
            if ret : 
                hap.append([gene_id, gene_name.replace('"', ''), int(fro), int(to), ret])
    if not hap : continue
    hap = sorted(hap, key = lambda x: int(x[2]))
    hap_str = hap[0][4]
    gene_arr = [hap[0][1] + ':' + str(hap[0][2]) + '..' + str(hap[0][3])]
    joints = []
    jointrg = []
    L = len(hap)
    if L == 1 :
        if ctg_strand == '-' : hap_str = rev_hap_str(hap_str)
        print(ctg, 'strand='+ctg_strand, 'vdj=' + vdj_tag, hap_str,
              '.', '-1..-1', hap[0][1], sep = '\t')
        continue

    for i in range(1, L) :
        gene_arr.append(hap[i][1] + ':' + str(hap[i][2]) + '..' + str(hap[i][3]))
        dis = hap[i][2] - hap[i-1][3] - 1
        sep = '~'
        c0 = dis < 25
        lab1 = hap[i-1][4]
        lab2 = hap[i][4]
        not_VV = hap[i-1][1][3] != 'V' or hap[i][1][3] != 'V'
        c1 = lab1[-1] != '>' and lab2[0] != '<' and dis < 100 and not_VV
        if c0 or c1 : 
            sep = '|' + str(dis) + '|'
            vdj_tag = 'T'
            if ctg_strand == '+' : joints.append(hap[i-1][1] + sep + hap[i][1])
            else : joints.append(hap[i][1] + sep + hap[i-1][1])
            jointrg.append(hap[i-1][2])
            jointrg.append(hap[i-1][3])
            jointrg.append(hap[i][2])
            jointrg.append(hap[i][3])
        hap_str += sep + lab2
    if not joints : joints = ['.']
    if ctg_strand == '-' : 
        hap_str = rev_hap_str(hap_str)
        gene_arr = gene_arr[::-1]
        joints = joints[::-1]
    if jointrg :
        joint_fro, joint_to = [min(jointrg), max(jointrg)]
    else :
        joint_fro = -1
        joint_to = -1
    print(ctg, 'strand='+ctg_strand, 'vdj=' + vdj_tag, hap_str,
          ','.join(joints), str(joint_fro) + '..' + str(joint_to), 
          ",".join(gene_arr), sep = '\t')
