import sys,re
import pprint as pp


Tbags = ['human', 'non-human', 
         'v', 'd', 'j', 'c', 
         'complete', 'partial']


blacklist = {}
blacklistfile = sys.argv[2]

with open(blacklistfile, "rt") as fp:
    while True:
        line = fp.readline()
        if not line: break
        allele,tmp = line.split('|')
        ctg,rg = tmp.split(':')
        rg0,rg1 = rg[:-1].split('-')
        if ctg not in blacklist : blacklist[ctg] = {}
        if allele not in blacklist[ctg] : blacklist[ctg][allele] = []
        blacklist[ctg][allele].append([int(rg0), int(rg1)])


inptags = []

if len(sys.argv) > 1 :
    inptags = sys.argv[1].split(',')
    for tag in inptags :
        if tag not in Tbags :
            exit('Error: ' + tag + ' is not supported, should be one of ' +
                 ','.join(Tbags) + '\n')


def rev_complement(seq) :
    rc_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    ret = []
    for x in seq[::-1] :
        if x in rc_map :
            ret.append(rc_map[x])
        else :
            ret.append(x)
    return(''.join(ret))


def read_id(line) :
    regex = re.compile(r'AC\s+(.*);')
    ret = regex.search(line)
    sid = ''
    if ret :
        sid = ret.group(1)
    return(sid)

def read_organism(line) :
    regex = re.compile(r'OS\s+(.*)\n')
    ret = regex.search(line)
    organism = ''
    if ret :
        organism = ret.group(1)
    return(organism)

def read_seq(buf) :
    llst = buf.split('\n')
    tag = False
    seq = []
    for line in llst :
        if tag :
            seq.append("".join(filter(lambda x: x.isalpha(), line)))
        if line[0:2] == 'SQ' :
            tag = True
    return("".join(seq))


def split_seq(seq, cut) :
    ret = ''
    i = 0
    while i < len(seq) :
        ret += seq[i:(i+cut)] + '\n'
        i += cut
    return(ret)

def is_partial(annot) :
    if 'partial' in annot : return(True)
    else : return (False)

def find_IMGT_allele(annot) :
    regex = re.compile(r'IMGT_allele="(.*?)"')
    ret = regex.search(annot)
    if ret :
        return(ret.group(1))
    else :
        return('')

def find_IMGT_gene(annot) :
    regex = re.compile(r'IMGT_gene="(.*?)"')
    ret = regex.search(annot)
    if ret :
        return(ret.group(1))
    else :
        return('')

def find_codon_start(annot) :
    regex = re.compile(r'codon_start=(\d)')
    ret = regex.search(annot)
    if ret :
        return(int(ret.group(1))-1)
    else :
        return('')

def find_gene_function(annot) :
    if '/functional' in annot: return('F')
    elif '/pseudo' in annot: return('P')
    elif '/ORF' in annot: return('ORF')
    else : return('U')


def check_arr_space(arrs) :
    s_arrs = sorted(arrs, key=lambda x: int(x[0]))
    space = 0
    a = s_arrs[0][0]
    for x in s_arrs :
        b = int(x[0])
        if b - a > space : space = b - a
        a = int(x[1])
    return(space)


def read_FT(buf, genetag):
    lines = buf.split('\n')

    complete_cut = 0
    if genetag == 'v' :
        key1s = ["V-GENE"]
        key2s = ["L-PART1", "L-PART2",
                 "V-HEPTAMER", "V-SPACER", "V-NONAMER",
                ]
        key3s = ["V-REGION"]
        complete_cut = 6
    elif genetag == 'd' :
        key1s = ["D-GENE"]
        key2s = ["5'D-NONAMER", "5'D-SPACER", "5'D-HEPTAMER",
                 "3'D-NONAMER", "3'D-SPACER", "3'D-HEPTAMER"]
        key3s = ["D-REGION"]
        complete_cut = 7
    elif genetag == 'j' :
        key1s = ["J-GENE"]
        key2s = ["J-NONAMER", "J-SPACER", "J-HEPTAMER"]
        key3s = ["J-REGION"]
        complete_cut = 4
    elif genetag == 'c' :
        key1s = ["C-GENE"]
        key2s = ["CH1","CH2","CH3","CH4","CHS", "CH5", "CH6", "CH7",
                 "CH1D","CH2D","CH2D2","CH2D3","CH3D","CH3D2","CH3D3",
                 "CH4D","CH4D2","CH4D3", "CH8", "CH9", "CH10","CHT","CHX",
                 "CL", "EX0", "EX4", "EX5", "EX6", "EX7", "EX8", "EX9",
                 "EX1","EX2","EX2A","EX2B","EX2R","EX2T","EX3","EX4UTR",
                 "H","I-EXON","M","M1","M2", "H1", "H2", "H3", "H4", "H5",
                 "STOP-CODON"]
        key2_exclude = ["ACCEPTOR-SPLICE", "INTRON", "TM-PART1",
                        "CONNECTING-REGIO", "DONOR-SPLICE", "C-GENE-UNIT",
                        "MISC_FEATURE", "TM-PART2", "CYTOPLASMIC-REGI",
                        "CL", "REPEAT_UNIT", "DELETION", "J-C-INTRON",
                        "INSERTION", "INT-DONOR-SPLICE", "TRANSMEMBRANE-RE",
                        "CO-PART2", "CO-PART3", "CO-PART4", "CO-PART1", 
                        "ENHANCER", "MUTATION", "GAP", "C-CLUSTER",
                        "PRIMER_BIND", "POLYA_SIGNAL", "3'UTR", "5'UTR",
                       ]
        key3s = ["C-REGION"]
    else :
        exit(genetag + " is not supported")

    L = len(lines)
    out0 = []
    gene_rgs = []
    key = ''
    annot = ''
    rc_tag = '+'
    regex = re.compile(r'<*(\d+)..>*(\d+)')
    for line in lines :
        if line[0:2] == "FT" :
            k_ret = line[3:21].split()
            if k_ret : 
                if key and rg:
                    out0.append([key, rg, rc_tag, annot])
                    annot = ''
                key = k_ret[0]
                if 'complement' in line : rc_tag = '-'
                else : rc_tag = '+'
                ret = regex.search(line)
                if ret :
                    rg = [int(ret.group(1)), int(ret.group(2))]
                    if key in key1s :
                        gene_rgs.append(rg)
                else :
                    rg = []
            else :
                annot += line[21:].strip()
    if key and rg:
        out0.append([key, rg, rc_tag, annot])


    out1 = []
    for rg in gene_rgs:
        tmp = []
        rec = {}
        partial = False
        gene_func = ''
        IMGT_allele = ''
        IMGT_gene = ''
        rc_tag = '+'
        new_rg = [99999999999, 0]
        for x in out0:
            if x[1][0] >= rg[0] and x[1][1] <= rg[1] :
                annot = x[3]
                if x[0] in key1s :
                    partial = is_partial(annot)
                    IMGT_allele = find_IMGT_allele(annot)
                    IMGT_gene = find_IMGT_gene(annot)
                    gene_func = find_gene_function(annot)
                    rc_tag = x[2]
                elif x[0] in key2s or x[0] in key3s:
                    if not partial : partial = is_partial(annot)
                    if not IMGT_allele: IMGT_allele = find_IMGT_allele(annot)
                    if not IMGT_gene : IMGT_gene = find_IMGT_gene(annot)
                    if not gene_func : gene_func = find_gene_function(annot)
                    if new_rg[0] > x[1][0] : new_rg[0] = x[1][0]
                    if new_rg[1] < x[1][1] : new_rg[1] = x[1][1]
                    tmp.append(x)

        rec['complete'] = not partial
        rec['gene_func'] = gene_func
        rec['IMGT_allele'] = IMGT_allele
        rec['IMGT_gene'] = IMGT_gene
        rec['rg'] = new_rg
        rec['rc_tag'] = rc_tag
        rec['elements'] = []
        tmp = sorted(tmp, key=lambda x:int(x[1][0]))
        for x in tmp:
            rec['elements'].append([x[0], x[1], find_codon_start(x[3])])
        # if the annotation does not contain 'partial' tag but it still miss
        # lots of elements, we will say it is a partial rather than a complete
        # gene

        if genetag != 'c' and rec['complete'] :
            if len(rec['elements']) != complete_cut :
                rec['complete'] = False
        
        if genetag == 'v' :
            arrs = []
            # no space should be between the V elements
            for x in rec['elements'] :
                if x[0] in ["L-PART2", "V-REGION",
                            "V-HEPTAMER", "V-SPACER","V-NONAMER"] :
                    arrs.append(x[1])
            if arrs and check_arr_space(arrs) > 1 :
                # only keep V region
                for x in rec['elements'] :
                    if x[0] == "V-REGION" :
                        rec['elements'] = [x]
                        rec['rg'] = x[1]
                        rec['complete'] = False
                        break
            #exclude extreme small V segments
            if rec['rg'][1] - rec['rg'][0] < 50 :
                rec = []
        if rec:
            out1.append(rec)

    return(out1)

def extract_data(buf, organism, Kbags) :
    meta_v = []
    meta_d = []
    meta_j = []
    meta_c = []
    meta = []
    seq = ''
    if 'human' in Kbags and 'human' not in organism :
        return([meta, seq])
    elif 'non-human' in Kbags and 'human' in organism :
        return([meta, seq])
    else :
        pass
    seq = read_seq(buf)

    if Kbags :
        if 'v' in Kbags : meta_v = read_FT(buf, 'v')
        if 'd' in Kbags : meta_d = read_FT(buf, 'd')
        if 'j' in Kbags : meta_j = read_FT(buf, 'j')
        if 'c' in Kbags : meta_c = read_FT(buf, 'c')
        meta = meta_v + meta_d + meta_j + meta_c

    if len(set(['v', 'd', 'j', 'c']).intersection(Kbags)) == 0 :
        meta_v = read_FT(buf, 'v')
        meta_d = read_FT(buf, 'd')
        meta_j = read_FT(buf, 'j')
        meta_c = read_FT(buf, 'c')
        meta = meta_v + meta_d + meta_j + meta_c

    if Kbags :
        meta0 = []
        if 'complete' in Kbags and 'partial' in Kbags :
            exit('Error: conflicted tags :' + 'partial' + ' VS ' + 'complete')
        elif 'complete' in Kbags :
            for x in meta :
                if x['complete'] : meta0.append(x)
            meta = meta0
        elif 'partial' in Kbags :
            for x in meta :
                if not x['complete'] : meta0.append(x)
            meta = meta0
        else :
            pass

    return([meta, seq])

buf = ''

for line in sys.stdin:
    if not line: break
    if line[0:2] == 'AC' :
        sid = read_id(line).replace(' ', '')
    if line[0:2] == 'OS' :
        organism = read_organism(line)
    buf += line
    if '//' in line :
        meta, seq = extract_data(buf, organism, inptags)

        buf = ''
        if not meta: continue

        for x in meta:
            if not sid or not x['IMGT_allele'] or \
               not x['rg'] or not x['rc_tag'] or not x['elements']:
                pass
            elif 'KIR' in x['IMGT_allele'] :
                pass
            else:
                total_length = x['rg'][1] - x['rg'][0] + 1
                rc_tag = x['rc_tag']
                allele_out = x['IMGT_allele'].replace(' ', '_')
                if sid in blacklist and allele_out in blacklist[sid] : continue
                allele_out += '|' + sid + ':' + str(x['rg'][0]) + '-' + str(x['rg'][1])

                ostr = '>' + allele_out
                ostr += ' ' + x['gene_func']
                fro = x['rg'][0]-1
                to = x['rg'][1]+1
                if rc_tag == '+' :
                    for y in x['elements'] :
                        ostr += ' ' + y[0] + '=' + str(y[1][0]-fro)
                        ostr += '..' + str(y[1][1]-fro)
                else :
                    x['elements'] = sorted(x['elements'], key = lambda y : y[1][1], reverse=True)
                    for y in x['elements'] :
                        ostr += ' ' + y[0] + '=' + str(to - y[1][1])
                        ostr += '..' + str(to - y[1][0])
                ostr += ' tot_length=' + str(total_length)
                ostr += '\n'

                oseq = seq[fro:(to-1)]
                if rc_tag == '-' :
                    ostr += split_seq(rev_complement(oseq), 80)
                else :
                    ostr += split_seq(oseq, 80)
                sys.stdout.write(ostr)
