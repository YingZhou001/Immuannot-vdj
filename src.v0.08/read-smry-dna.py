import sys,re
import pprint as pp


def is_constant_gene(gene) :
    c_genes = ["IGHA","IGHA1","IGHA2","IGHA4","IGHA5","IGHD","IGHD1","IGHDD","IGHDD1P","IGHDD2P","IGHDD3P","IGHE","IGHEP1","IGHG","IGHG1","IGHG1A","IGHG1B","IGHG2","IGHG2A","IGHG2B","IGHG2C","IGHG3","IGHG3A","IGHG3B","IGHG3C","IGHG4","IGHG4A","IGHG5","IGHG5-1","IGHG5-2","IGHG6","IGHG6-1","IGHG6-2","IGHG7","IGHGP","IGHM","IGHM1","IGHM2","IGHM6","IGHMD","IGHO","IGHT1","IGHT1D","IGHT2","IGHT2D","IGHT4","IGHT5","IGHZ","IGIC1","IGIC2","IGIC3","IGIC4","IGIC5","IGKC","IGKC1","IGKC2","IGLC","IGLC1","IGLC10","IGLC11","IGLC2","IGLC2A","IGLC3","IGLC4","IGLC4P","IGLC5","IGLC6","IGLC7","IGLC8","IGLC9","IGLC/OR17-1","IGLC/OR17-2","IGLC/OR17-3","IGLC/OR22-1","IGLC/OR22-2","TRAC","TRBC1","TRBC2","TRBC3","TRDC","TRGC","TRGC1","TRGC2","TRGC2S1","TRGC3","TRGC4","TRGC5","TRGC6","TRGC7","TRGC8"]
    ret = gene.replace('"', '') in c_genes
    return(ret)

db = {}
read_info = {}


for line in sys.stdin:
    if not line : break
    if line[0] == '#': continue
    llst = line.split('\t')
    feature = llst[2]
    if feature in ['transcript', 'transcript(partial)'] : continue
    attr = {}
    for x in llst[-1][:-1].split('; ') :
        k,v = x.split(' ')
        attr[k] = v
    ctg = llst[0]
    fro = int(llst[3])
    to = int(llst[4])
    strand = llst[6]
    size = to - fro + 1

    if ctg not in db : db[ctg] = []
    gene_name = attr['gene_name'].replace('"', '')
    if is_constant_gene(gene_name) :
        if feature not in ['gene', 'gene(partial)'] : continue
    
    if feature in ['gene', 'gene(partial)'] : 
        db[ctg].append([fro, to, strand, gene_name, []])
    else :
        db[ctg][-1][4].append([fro,to,feature])


def add_symbol(da) :
    for i in range(len(da)) :
        strand = da[i][2]
        gene = da[i][3]
        x = da[i][4]
        sstr = []
        if not is_constant_gene(gene) :
            for y in x :
                feature_tag = convert_feature(y[2], strand, gene)
                sstr.append(feature_tag)
        else :
            feature_tag = convert_feature('', strand, gene)
            sstr.append(feature_tag)
        sstr = ''.join(sstr)
        sstr = sstr.replace('SSR', '<')
        sstr = sstr.replace('SR', '<')
        sstr = sstr.replace('RSS', '>')
        sstr = sstr.replace('RS', '>')
        sstr = sstr.replace('R', '') # remove independent heptamer
        if 'V' in sstr.upper() : pass
        elif 'D' in sstr.upper() : pass
        elif 'J' in sstr.upper() : pass 
        elif is_constant_gene(gene) : pass
        else :
            da[i][3] = '.'
        da[i][4] = ''.join(sstr)
    return(da)


def rm_standby_L1(da) :
    if not da : return(da)
    L = len(da)
    idx_to_rm = []
    for i in range(L) :
        if len(da[i][4]) == 1 and 'L-PART1' in da[i][4][0] :
            idx_to_rm.append(i)
    newda = []
    for i in range(L) :
        if i not in idx_to_rm :
            newda.append(da[i])
    return(newda)



def check_strand_and_region(da) :
    pos_l = 0
    neg_l = 0
    reg_l = {}
    for fro, to, strand, gene_name, tmp in da :
        length = int(to) - int(fro)
        reg = gene_name[0:3]
        if reg not in reg_l : reg_l[reg] = 0
        reg_l[reg] += length

        if strand == '+' : pos_l += length
        else : neg_l += length
    ret_reg = ''
    ret_reg_l = 0
    for reg in reg_l :
        if ret_reg_l < reg_l[reg] : 
            ret_reg = reg
            ret_reg_l = reg_l[reg]
    if pos_l >= neg_l : ret_strand = '+'
    else : ret_strand = '-'
    return([ret_strand, ret_reg])

def check_constant_gene(da) :
    if not da : return('')
    tmp = {}
    for fro,to,gene in da :
        if gene not in tmp :
            tmp[gene] = 0
        tmp[gene] += to - fro + 1
    ret_gene = ''
    ret_l = 0
    for gene in tmp :
        if tmp[gene] > ret_l :
            ret_gene = gene
            ret_l = tmp[gene]
    return(ret_gene)

def rev_hap_str(hap_str) :
    return()
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


def convert_feature(feature, strand, gene):
    sym_map = {'L-PART1': '^', 'L-PART2': '^',
               'V-REGION': 'V', 'D-REGION': 'D', 'J-REGION': 'J',
               'J-NONAMER':'S', 'J-SPACER':'S', 'J-HEPTAMER': 'R',
               'V-NONAMER':'S', 'V-SPACER':'S', 'V-HEPTAMER': 'R',
               '5\'D-NONAMER':'S', '5\'D-SPACER':'S', '5\'D-HEPTAMER': 'R',
               '3\'D-NONAMER':'S', '3\'D-SPACER':'S', '3\'D-HEPTAMER': 'R', 
              }
    if feature in sym_map :
        ret = sym_map[feature]
    elif is_constant_gene(gene):
        ret = gene[3:]
    else :
        ret = '.'
    if strand == '-' and ret not in ['R', 'S', '^', '$'] :
        ret = ret.lower()
    return(ret)

strad_rev_map = {'-':'+', '+':'-'}

hap_str = {}

for ctg in db:
    vdj_tag = 'F'
    [ctg_strand, ctg_reg] = check_strand_and_region(db[ctg])
    if ctg_strand == '-' :
        db[ctg] = sorted(db[ctg], key = lambda x : x[0], reverse=True)
        for i in range(len(db[ctg])) :
            db[ctg][i][2] = strad_rev_map[db[ctg][i][2]]
            db[ctg][i][4] = sorted(db[ctg][i][4], key = lambda x : x[0], reverse=True)
    else :
        db[ctg] = sorted(db[ctg], key = lambda x : x[0], reverse=False)

    #db[ctg] = rm_standby_heptamer(db[ctg])
    db[ctg] = rm_standby_L1(db[ctg])
    db[ctg] = add_symbol(db[ctg])

    if not db[ctg] : continue

    if ctg not in hap_str : hap_str[ctg] = []

    da = db[ctg]
    L = len(da)
    all_genes = []
    old_fro, old_to, strand, old_gene, old_tag = da[0]
    hap_str[ctg].append(old_tag)
    all_genes.append(old_gene)

    vj_joint = []
    dj_joint = []
    dd_joint = []
    vd_joint = []
    jj_joint = []
    
    for i in range(1,L) :
        fro,to,strand,gene,tag = da[i]
        dis = min(abs(fro - old_to), abs(old_fro - to))
        g1 = old_gene
        g2 = gene
        f1 = old_tag.upper()
        f2 = tag.upper()
        rg0 = min([old_fro, old_to, fro, to])
        rg1 = max([old_fro, old_to, fro, to])
        rg = str(rg0) + '..' + str(rg1)

        is_const = is_constant_gene(g1) or is_constant_gene(g2)

        no_rss = '>' not in f1 and '<' not in f2
        is_vj = ('V' in f1 and 'J' in f2) and '<' not in f2
        is_vj = is_vj or ('J' in f1 and 'V' in f2) and '>' not in f1
        is_dd = ('D' in f1.upper() and 'D' in f2.upper()) and  no_rss
        is_dd = is_dd and not is_const
        is_dj = 'D' in f1 and 'J' in f2 and no_rss
        is_dj = is_dj or ( 'J' in f1 and 'D' in f2 and no_rss )
        is_dj = is_dj and not is_const
        is_vd = 'V' in f1 and 'D' in f2 and '<' not in f2
        is_vd = is_vd or 'D' in f1 and 'V' in f2 and '>' not in f1
        is_vd = is_vd and not is_const
        is_jj = 'J' in f1.upper() and 'J' in f2.upper() and \
                not ('>' in f1 or '<' in f2 or '<' in f1 or '>' in f2)
        is_rss_joint = '<' in f1 and '>' in f2
        if is_vj :
            if ctg_reg != 'IGH' :
                if dis >= 50 : 
                    hap_str[ctg].append('~')
                    all_genes.append('~')
                else :
                    hap_str[ctg].append('|')
                    all_genes.append('|')
                    vj_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
            else :
                if dis >= 300 :
                    hap_str[ctg].append('~')
                    all_genes.append('~')
                else :
                    hap_str[ctg].append('|')
                    all_genes.append('|')
                    vj_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
        elif is_dd :
            if dis >= 50 :
                hap_str[ctg].append('~')
                all_genes.append('~')
            else :
                hap_str[ctg].append('|')
                all_genes.append('|')
                dd_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
        elif is_dj :
            if dis >= 50 :
                hap_str[ctg].append('~')
                all_genes.append('~')
            else :
                hap_str[ctg].append('|')
                all_genes.append('|')
                dj_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
        elif is_vd :
            if dis >= 50 :
                hap_str[ctg].append('~')
                all_genes.append('~')
            else :
                hap_str[ctg].append('|')
                all_genes.append('|')
                vd_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
        elif is_jj :
            if dis >= 50 :
                hap_str[ctg].append('~')
                all_genes.append('~')
            else :
                hap_str[ctg].append('|')
                all_genes.append('|')
                jj_joint.append(g1 + '|' + str(dis) + '|' + g2 + ',' + rg)
        elif is_rss_joint :
            if dis >= 10 :
                hap_str[ctg].append('~')
                all_genes.append('~')
            else :
                hap_str[ctg].append('|')
                all_genes.append('|')
        elif dis < 10 :
            hap_str[ctg].append('|')
            all_genes.append('|')
        else :
            hap_str[ctg].append('~')
            all_genes.append('~')

        hap_str[ctg].append(tag)
        all_genes.append(gene)

        old_to = to
        old_fro = fro
        old_tag = tag
        old_gene = gene
        
    if not vj_joint:  vj_joint.append('.')
    if not vd_joint:  vd_joint.append('.')
    if not dd_joint:  dd_joint.append('.')
    if not dj_joint:  dj_joint.append('.')
    if not jj_joint:  jj_joint.append('.')

    out_hapstr = ''.join(hap_str[ctg])
    out_genestr = ''.join(all_genes)

    if len(all_genes) > 0 :
        print(ctg,ctg_reg, out_genestr, out_hapstr,
              'vj=' + ';'.join(vj_joint), 
              'vd=' + ';'.join(vd_joint), 
              'dd=' + ';'.join(dd_joint), 
              'dj=' + ';'.join(dj_joint), 
              'jj=' + ';'.join(jj_joint), 
              sep='\t')
