import sys,gzip,re


# load the meta info for each gene

gene_db = {}

with gzip.open(sys.argv[1], 'rt') as fp:
    while True:
        line = fp.readline()
        if not line: break
        if line[0] != '>' : continue
        llst = line[1:].split()
        source = llst[0]
        gene_db[source] = []
        for x in llst:
            if '..' in x :
                ele,fro,to = re.split('=|\.\.', x)
                # change 1-based into 0-based
                gene_db[source].append([ele, int(fro), int(to)])


def overlap(a,b, c,d) :
    # determine whether [c, d) being full/partial/no included [a,b)
    if b <= c or a >= d :
        return('')
    elif a <= c and d <= b :
        return('f')
    else :
        return('p')

def parsing_ds(ds) :
    matches = re.findall(r'([:\+\-\*]{1})([\[\]\da-z]+)', ds)
    out = [[m[0], m[1]] for m in matches]
    return(out)

def cal_ds_nm(ds) :
    mat = parsing_ds(ds)
    ds_nm = 0
    for m in mat :
        tag = m[0]
        info = m[1]
        if tag == '*' : ds_nm += 1
        elif tag in ['+', '-'] :
            if '[' in info and len(info) == 3 :
                continue
            else :
                ds_nm += 1
    return(ds_nm)

def extract_ds(paf_rec) :
    regex = re.compile(r'ds:Z:(\S*)')
    ret = regex.search(paf_rec)
    if ret :
        ds_str = ret.group(1)
    else :
        ds_str = ''
    return(ds_str)

def query_to_ref_map(qx, query_rg, ref_rg, strand, cs) :
    # qx, 1-based coordinate
    # query_rg, 0-based
    # ref_rg, 0-based
    # process in 1-based coordinate
    # output, 1 based
    if '[' in cs : exit("please use cs string as input")
    if strand == '+' :
        if qx <= query_rg[0] : return(ref_rg[0]+1)
        elif qx >= query_rg[1] : return(ref_rg[1])
        qx = qx - query_rg[0]
    elif strand == '-' :
        if qx <= query_rg[0] : return(ref_rg[1])
        elif qx >= query_rg[1] : return(ref_rg[0]+1)
        qx = query_rg[1] - qx + 1
    else :
        exit('no valid strand information')

    mat = parsing_ds(cs)
    old_ref_p = 0
    old_qry_p = 0
    qy = -1
    for tag,op in mat :
        if tag == ':' :
            ref_inc = int(op)
            qry_inc = int(op)
        elif tag == '*' :
            ref_inc = 1
            qry_inc = 1
        elif tag == '+' :
            qry_inc = len(op)
            ref_inc = 0
        elif tag == '-' :
            ref_inc = len(op)
            qry_inc = 0
        cur_ref_p = old_ref_p + ref_inc
        cur_qry_p = old_qry_p + qry_inc

        dif = qx - old_qry_p

        if dif <= qry_inc :
            if tag in [':', '*'] :
                qy = old_ref_p + dif
            elif tag == '+' :
                qy = old_ref_p
            #qy = old_ref_p + dif
            break
        old_ref_p = cur_ref_p
        old_qry_p = cur_qry_p 

    if qy < 0: 
        print("error: in add-elements.py", qx, query_rg, ref_rg, strand, cs, file=sys.stderr)
        exit()
    qy += ref_rg[0]
    return(qy)



print("##gff-version 2")


gff_tag = "Immuannot-VDJ"
idx = 1
for line in sys.stdin:
    if not line : break
    llst = line.split()
    ctg = llst[0]
    fro = int(llst[1])
    to = int(llst[2])
    gene = llst[3].split('*')[0]
    strand = gene[0]
    gene = gene[1:]
    qual = round(float(llst[5].split(',')[0][:-1].split('=')[1]), 2)
    paf = line[:-1].split('@')[1]
    ds_str = extract_ds(paf)

    if ds_str :
        ds_nm = cal_ds_nm(ds_str)
        cs_str = re.sub('\[|\]', '', ds_str)
    else :
        ds_nm = -9
        cs_str = ''
    llst2 = paf.split()
    source = llst2[0]
    allele = source.split('|')[0].split('=')[0]
    map_fro = int(llst2[2])
    map_to = int(llst2[3])
    if qual < 100 : 
        gene_field_tag = 'gene(partial)'
    else : 
        gene_field_tag = 'gene'
    gene_id = 'GNM' + str(idx).rjust(5, '0')
    idx += 1
    gene_row = [ctg, gff_tag, gene_field_tag, str(fro+1), str(to), '.', strand, '.']
    attr = 'gene_id "' + gene_id + '";'
    attr += ' gene_name "' + gene + '";'
    attr += ' template_source "' + source + '";'
    attr += ' mapping_rate ' + str(qual/100) + ';'
    attr += ' ds_nm ' + str(ds_nm)
    gene_row.append(attr)
    print('\t'.join(gene_row))

    #check elements
    eles = []
    if source in gene_db :
        query_rg = [map_fro, map_to]
        ref_rg = [fro, to]
        cs = cs_str
        #if gene == 'IGHG3':  print(cs, file=sys.stderr)
        for x in gene_db[source] :
            ele, ele_fro, ele_to = x
            #if gene == 'IGHG3':  print(x, query_rg, ref_rg, strand, file=sys.stderr)
            read_fro = query_to_ref_map(ele_fro, query_rg, ref_rg, strand, cs)
            read_to = query_to_ref_map(ele_to, query_rg, ref_rg, strand, cs)
            read_fro, read_to = sorted([read_fro, read_to])
            prop = (read_to + 1 - read_fro)/(ele_to - ele_fro + 1)
            if prop < 0.05 or read_fro == read_to : continue
            ele_row = [ctg, gff_tag, ele, str(read_fro), str(read_to), '.', strand, '.']
            attr = 'gene_id "' + gene_id + '";'
            attr += ' gene_name "' + gene + '";'
            attr += ' length_ratio ' + str(round(prop, 2))
            ele_row.append(attr)
            eles.append(ele_row)
    eles = sorted(eles, key = lambda x: int(x[3]))
    for ele in eles :
        print('\t'.join(ele))
