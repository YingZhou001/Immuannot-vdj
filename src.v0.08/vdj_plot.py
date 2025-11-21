import sys,re,gzip
import pprint as pp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser(description='An lightwieght program to massively vissulize VDJ recombination and more')

parser.add_argument('--paf', type=str, 
                    help='(required) alignment file from minimap2')
parser.add_argument('--refbed', type=str, 
                    help='(required) gene annotation for the reference')
parser.add_argument('--tarbed', type=str, default="",
                    help='gene annotation for the target')
parser.add_argument('--refmeta', type=str, default="",
                    help='information to add to the reference panel')
parser.add_argument('--autorev', action='store_true',
                   help='automatically reverse the strand')
parser.add_argument('--order', action='store_true',
                   help='read plot in the order of meta file')
parser.add_argument('--tarmeta', type=str, default="",
                    help='information to add to the target panel')
parser.add_argument('--outpref', type=str, default="out-vdj",
                    help='output prefix')
parser.add_argument('--layout', type=str, default="5,2",
                    help='layout of each page in nrow,ncol')
parser.add_argument('--size', type=str, default="25,15",
                    help='size of each page in width,height')
parser.add_argument('--gapcompress', type=str, default="1000",
                    help='gap compression ratio')
parser.add_argument('--gapfill', type=str, default="10",
                    help='gap fill cutoff')



args = parser.parse_args()

paffile = args.paf
refbedfile = args.refbed
tarbedfile = args.tarbed
refmetafile = args.refmeta
tarmetafile = args.tarmeta
autorev_tag = args.autorev
order_tag = args.order
gapratio = int(args.gapcompress)
layout = [int(x) for x in args.layout.split(',')]
figsize = [int(x) for x in args.size.split(',')]
gapfill = int(args.gapfill)
outpref = args.outpref

# process paf
def parse_cigar(cigar) :
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar)]

def merge_gap(splt_mat, gap_fill=gapfill) :
    new_mat = []
    old_ref_fro = -gap_fill
    old_ref_to = -gap_fill
    old_query_fro = -gap_fill
    old_query_to = -gap_fill

    for query_fro, query_to, ref_fro, ref_to, op in splt_mat :
        ref_l = ref_to - ref_fro
        query_l = query_to - query_fro

        if op in ['D' , 'I' , 'N']  : continue

        ref_dis = abs(ref_fro + ref_to - old_ref_fro - old_ref_to)/2
        ref_dis -= (ref_to - ref_fro + old_ref_to - old_ref_fro)/2 

        query_dis = abs(query_fro + query_to - old_query_fro - old_query_to)/2
        query_dis -= (query_to - query_fro + old_query_to - old_query_fro)/2 

        if ref_dis < gap_fill and query_dis < gap_fill :
            new_mat[-1][0] = min(query_fro, old_query_fro)
            new_mat[-1][1] = max(query_to, old_query_to)
            new_mat[-1][2] = min(ref_fro, old_ref_fro)
            new_mat[-1][3] = max(ref_to, old_ref_to)
        else :
            new_mat.append([query_fro, query_to, ref_fro, ref_to])
        
        old_query_fro,old_query_to,old_ref_fro,old_ref_to = new_mat[-1][0:4]
    return(new_mat)
    


def b_split(query_rg, ref_rg, cigar, strand) :
    # return a matrix and each row is
    # query_block, ref_block, strand
    cg_mat = []
    cur_query_pos = 0
    cur_ref_pos = 0
    for l,op in parse_cigar(cigar) :
        cg_mat.append([cur_query_pos, cur_ref_pos, l, op])
        if op in ['M', '=', 'X', 'D', 'N'] :
            cur_ref_pos = cur_ref_pos+l
        if op in ['M', '=', 'X', 'I'] :
            cur_query_pos = cur_query_pos+l
        if op in ['S', 'H', 'P'] :
            exit("not support 'SHP' in cigar!")

    out = []
    ref_shift = ref_rg[0]
    query_shift = 0

    for q,r,l,op in cg_mat :
        if op in ['M', '=', 'X'] :
            ref_fro = r + ref_shift
            ref_to = ref_fro + l
            query_fro = q + query_shift
            query_to = query_fro + l
        elif op in ['D', 'N']:
            ref_fro = r + ref_shift
            ref_to = ref_fro + l
            query_fro = q + query_shift
            query_to = query_fro
        elif op == 'I' :
            ref_fro = r + ref_shift
            ref_to = ref_to
            query_fro = q + query_shift
            query_to = query_fro + l

        out.append([query_fro, query_to, ref_fro, ref_to, op])
        #print([query_fro, query_to, ref_fro, ref_to, l, op])

    L = len(out)
    if strand == '+' :
        for i in range(L) :
            out[i][0] = out[i][0] + query_rg[0]
            out[i][1] = out[i][1] + query_rg[0]
    elif strand == '-' :
        for i in range(L) :
            a,b = out[i][0:2]
            out[i][0] = query_rg[1] - b
            out[i][1] = query_rg[1] - a
    return(out)

def is_gz_file(filepath):
    with open(filepath, 'rb') as fp:
        return fp.read(2) == b'\x1f\x8b'

def load_bed(bedfile) :
    bed = {}
    ctg_list = []
    if not bedfile :
        return([bed, ctg_list])
    gz_tag = is_gz_file(bedfile)
    if gz_tag :
        fp = gzip.open(bedfile, 'rt')
    else :
        fp = open(bedfile, 'rt')
    while True:
        line = fp.readline()
        if not line: break
        llst = line.split()
        ctg = llst[0]
        fro = int(llst[1])
        to = int(llst[2])
        ele = llst[3]
        gene = llst[4]
        gid = llst[5]
        if ctg not in bed : bed[ctg] = []
        bed[ctg].append([fro, to, ele, gene, gid])
        if ctg not in ctg_list :
            ctg_list.append(ctg)
    fp.close()
    return([bed, ctg_list])

def load_meta(metafile) :
    meta = {}
    if not metafile : 
        return meta
    gz_tag = is_gz_file(metafile)
    if gz_tag :
        fp = gzip.open(metafile, 'rt')
    else :
        fp = open(metafile, 'rt')
    while True:
        line = fp.readline()
        if not line: break
        llst = line.split()
        ctg = llst[0]
        info = ','.join(llst[1:])
        if ctg not in meta : meta[ctg] = []
        meta[ctg].append(info)
    fp.close()
    return(meta)
   
def load_meta_order(metafile) :
    meta_order = []
    if not metafile : 
        return meta_order
    gz_tag = is_gz_file(metafile)
    if gz_tag :
        fp = gzip.open(metafile, 'rt')
    else :
        fp = open(metafile, 'rt')
    while True:
        line = fp.readline()
        if not line: break
        llst = line.split()
        ctg = llst[0]
        meta_order.append(ctg)
    fp.close()
    return(meta_order)
 


def load_paf(paffile) :
    paf = {}
    keys = []
    cg_regex = re.compile(r'cg:Z:(.*?)\s')
    gz_tag = is_gz_file(paffile)
    if gz_tag :
        fp = gzip.open(paffile, 'rt')
    else :
        fp = open(paffile, 'rt')
    while True:
        line = fp.readline()
        if not line: break
        llst = line.split()
        qctg, qlen, qfro, qto, strand, ref, rlen, rfro, rto =llst[0:9]
        if not tar_ctg_list :
            if ref not in ref_ctg_list : continue
        elif not (qctg in tar_ctg_list and ref in ref_ctg_list) :
            continue
        key = qctg + ' ' + ref
        if key.count(' ') > 1 : exit('Format error: " "(blank) sould not in the contig name')
        if key not in paf : 
            paf[key] = {'length': [int(qlen), int(rlen)], 'mapping':[], 'rev_tag':[0,0]}
            keys.append(key)
        ret = cg_regex.search(line)
        if ret:
            cigar = ret.group(1)
        else :
            exit("no cigar string found, please re-run the alignment")
        query_rg = [int(qfro), int(qto)]
        ref_rg = [int(rfro), int(rto)]
        split_mat = b_split(query_rg, ref_rg, cigar, strand)
        mat = merge_gap(split_mat)

        for arr in mat :
            qfro1,qto1, rfro1, rto1 = arr
            if abs(qfro1-qto1) > 5 and abs(rfro1 - rto1) > 5:
                paf[key]['mapping'].append([qfro1, qto1, strand, rfro1, rto1])
                if strand == '+' : paf[key]['rev_tag'][0] += qto1 - qfro1
                elif strand == '-' : paf[key]['rev_tag'][1] += qto1 - qfro1
    fp.close()
    return([paf, keys])

def build_compression_map(tab_arr, flank = 50, ratio = 1000) :
    # [[fro,to],[fro,ro], ...]
    # merge intervals within flank=50bp
    tmp_sorted = sorted(tab_arr, key = lambda x: x[0])
    new_tab = [tmp_sorted[0]]
    for a,b in tmp_sorted[1:] :
        a0 = new_tab[-1][0]
        b0 = new_tab[-1][1]
        if a <= b0 + flank:
            new_tab[-1][1] = max(b, new_tab[-1][1])
        else :
            new_tab.append([a,b])
    bc_map = {'x':[0], 'y':[0], 'r':[1/ratio]}
    for a,b in new_tab :
        x0 = bc_map['x'][-1]
        y0 = bc_map['y'][-1]
        r0 = bc_map['r'][-1]
        x1 = max(a - flank, 1)
        y1 = y0 +( x1 - x0) * r0
        r1 = 1
        x2 = b + flank
        y2 = y1 + (x2 - x1) * r1
        r2 = 1/ratio
        bc_map['x'].append(x1)
        bc_map['y'].append(y1)
        bc_map['r'].append(r1)
        bc_map['x'].append(x2)
        bc_map['y'].append(y2)
        bc_map['r'].append(r2)
    
    return(bc_map)

def bi_int_search(arr, x):
    if x <= min(arr) : return(0)
    elif x >= max(arr) : return(len(arr)-1)
    else : pass

    # Define the search bounds
    left = 0
    right = len(arr) - 1  

    while left < right - 1 :
        mid = (left + right ) // 2  
        if arr[mid] <= x : 
            left = mid 
        elif arr[mid] > x : 
            right = mid  
    return(left)


def coordinate_compress(bc_map, pos_arr) :
    out_arr = []
    for x in pos_arr :
        idx = bi_int_search(bc_map['x'], x)
        x0 = bc_map['x'][idx]
        y0 = bc_map['y'][idx]
        r0 = bc_map['r'][idx]
        y = int(y0 + (x-x0) * r0)
        out_arr.append(y)
    return(out_arr)

def check_overlap(x, xs) :
    # check interval x is overlapped with any interval in xs
    x0,x1 = x
    for y0,y1 in xs :
        if y0 > x1 or y1 < x0 : 
            continue
        else :
            return(True)
    return(False)

# plot functions
def add_ele_box (ax, rg, h0, h1, fc = 'blue', ec ='black') :
    x = rg[0]
    y = h0
    rect_width = rg[1] - rg[0]
    rect_height = h1-h0
    rect = patches.Rectangle((x, y), rect_width, rect_height, fc=fc, ec=ec)
    ax.add_patch(rect)

def add_mid_line(ax, rg, h0, h1, lw = 2, color = "black") :
    x = rg
    y = [(h0 + h1) / 2] *2
    ax.plot(x, y, lw=lw, color = color)

def impute(x, p1, p2) :
    # known x and impute y on the line between point p1 and p2
    x1,y1 = p1
    x2,y2 = p2
    r = (y1 - y2) / (x1 - x2)
    y = y1 + (x - x1) * r
    return(y)

def add_mapping(ax, rg0, rg1, h0, h1, color='lightblue', alpha=0.5) :
    x = [rg0[0], rg0[1], rg1[1], rg1[0], rg0[0]]
    y = [h0, h0, h1, h1, h0]
    ax.fill(x, y, color=color, alpha=alpha)

def add_text(ax, text, x, y, size, ha='center', va='center') :
    ax.text(x, y, text, fontsize=size, color='black',
             ha=ha, va=va)
             #bbox=dict(facecolor='white', alpha=0.5))

def get_h(h = [0, 1], margin = [1, 1, 3, 1, 1], gap_ratio = 0.05):
    ml = len(margin)
    totl = sum(margin)
    toth = h[1] - h[0]
    out = []
    cur_h = h[0]
    for x in margin:
        h0 = cur_h
        inc = x/totl * toth
        cur_h = h0 + inc
        h0 = h0 + inc * gap_ratio/2
        h1 = cur_h - inc * gap_ratio/2
        out.append([h0, h1])
    return(out)

def ele_col(ele) :
    cols = ["gray", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"]
    col_map = {
        "5'D-NONAMER": cols[1],
        "5'D-SPACER": cols[2],
        "5'D-HEPTAMER":  cols[3],
        "D-REGION":  cols[4],
        "3'D-HEPTAMER": cols[5],
        "3'D-SPACER": cols[6],
        "3'D-NONAMER": cols[7],
        "J-NONAMER": cols[1],
        "J-SPACER": cols[2],
        "J-HEPTAMER": cols[3],
        "J-REGION": cols[4],
        "L-PART1": cols[8],
        "L-PART2": cols[9],
        "V-REGION": cols[4],
        "V-HEPTAMER": cols[5],
        "V-SPACER": cols[6],
        "V-NONAMER": cols[7],
        "CH1" : cols[4],
        "CH2" : cols[4],
        "CH3" : cols[4],
        "CH3-CHS" : cols[4],
        "CH4" : cols[4],
        "CH4-CHS" : cols[4],
        "CHS" : cols[4],
        "C-REGION" : cols[4],
        "D-REGION" : cols[4],
        "EX1" : cols[4],
        "EX2" : cols[4],
        "EX3" : cols[4],
        "EX4UTR" : cols[4],
        "H" : cols[4],
        "H-CH2" : cols[4],
        "M" : cols[4],
        "M1" : cols[4],
        "M2 " : cols[4]
    }
    if ele in col_map: return(col_map[ele])
    else : return(cols[0])

def plot_one_pair(da, ax, ctg, ref, sub_tarbed, sub_refbed,
                  tar_info = '', ref_info = '',
                  h = [0,1], gapratio = 1000,
                  margin = [0.5, 1, 0.5, 3, 0.5, 1, 0.5],
                  auto_strand_rev = autorev_tag,
                 ):
    flip = {'+':'-', '-':'+'}
    ctg_len = da['length'][0]
    ctg_arr = []
    ref_arr = []
    ref_rg = [99999999999, 0]
    if auto_strand_rev :
        rev_tag = da['rev_tag'][1] > da['rev_tag'][0] * 2
    else :
        rev_tag = False
    #pp.pprint(da['mapping'])
    da_L = len(da['mapping'])
    for i in range(da_L) :
        qfro,qto,strand,rfro,rto = da['mapping'][i]
        if rev_tag :
            qfro,qto = ctg_len - qto, ctg_len - qfro
            da['mapping'][i][0] = qfro
            da['mapping'][i][1] = qto
            da['mapping'][i][2] = flip[strand]
        ctg_arr.append([qfro, qto])
        ref_arr.append([rfro, rto])
        ref_rg[0] = min(ref_rg[0], rfro)
        ref_rg[1] = max(ref_rg[1], rto)

    #print(sorted(ctg_arr))
    ctg_map = build_compression_map(ctg_arr, ratio = gapratio)
    #pp.pprint(ctg_map)
    ref_map = build_compression_map(ref_arr, ratio = gapratio)
    new_da = []
    refshift = 999999999
    for qfro,qto,strand,rfro,rto in da['mapping'] :
        #qfro,qto = coordinate_compress(ctg_map, [qfro,qto])
        rfro,rto = coordinate_compress(ref_map, [rfro,rto])
        new_da.append([qfro,qto,rfro,rto,strand])
        if refshift > rfro: refshift = rfro

    if sub_tarbed :
        for i in range(len(sub_tarbed)) :
            qfro, qto = sub_tarbed[i][0:2]
            if rev_tag :
                qfro,qto = ctg_len - qto, ctg_len - qfro
            #qfro,qto = coordinate_compress(ctg_map, [qfro,qto])
            sub_tarbed[i][0] = qfro
            sub_tarbed[i][1] = qto
    #new_ctg_len = coordinate_compress(ctg_map, [ctg_len])[0]
    new_ctg_len = ctg_len
    new_sub_refbed = []
    for i in range(len(sub_refbed)) :
        qfro, qto, ele, gene, gid = sub_refbed[i]
        if qto < ref_rg[0] -50  or qfro > ref_rg[1] + 50 : continue
        qfro, qto = coordinate_compress(ref_map, [qfro, qto])
        if abs(qfro -qto) < 3 : continue
        new_sub_refbed.append([qfro, qto, ele, gene, gid])

    
    height = h[1] - h[0]
    ax.set_ylim(h[0] - height*0.05, h[1] + height*0.05)
    h_cnm, h_ctext, h_cbox, h_map, h_rbox, h_rtext, h_rnm = get_h(h, margin)
    # add ref and target ctg name
    size = 10
    x = 0
    y = h_cnm[1]
    add_text(ax, ' '.join([ctg, str(ctg_len) + 'bp', tar_info]),x, y, size, ha='left', va='top')
    y = h_rnm[0]
    add_text(ax, ' '.join([ref, ref_info]), x, y, size, ha='left', va='bottom')

    #print("new mapping")
    ref_intervals = []
    for x in new_da:
        rg0 = x[0:2]
        rg1 = [y - refshift for y in x[2:4]]
        ref_intervals.append(rg1)
        h0,h1 = h_map
        strand = x[4]
        color=''
        if strand == '-' : 
            color='grey'
            rg1 = rg1[::-1]
        else : 
            color='blue'
        #print(rg0, rg1)
        add_mapping(ax, rg0, rg1, h0, h1, color=color)
    #print("tar annotation")
    tar_genes = {}
    rg = [0, new_ctg_len]
    h0,h1 = h_cbox
    #print(h_cbox)
    # plot bar for the reads
    add_ele_box (ax, rg, h0, h1, fc='lightblue', ec='lightblue')
    if sub_tarbed: 
        for x in sub_tarbed:  
            #print(x)
            rg = x[0:2]
            ele = x[2]
            col = ele_col(ele)
            add_ele_box (ax, rg, h0, h1, fc=col, ec=col)
            gene_k = x[3] + ' ' + x[4]
            if gene_k not in tar_genes : tar_genes[gene_k] = []
            tar_genes[gene_k].append(rg[0])
            tar_genes[gene_k].append(rg[1])
    
    h0,h1 = h_ctext
    textrow_i = 0
    textrow_inc = (h1 - h0)/4
    text_size = 10
    ha = 'center'
    for gene_k in tar_genes:
        pos_arr = tar_genes[gene_k]
        rg = [min(pos_arr), max(pos_arr)]
        x = sum(rg)/2
        if textrow_i % 2 == 0 :
            y = h0 + textrow_inc
        else :
            y = h0 + textrow_inc * 3
        text = gene_k[3:].split()[0]

        add_text(ax, text, x, y, size, ha=ha)
        textrow_i += 1
        add_mid_line(ax, rg, h_cbox[0], h_cbox[1], lw=2)

    #print("ref annotation")
    ref_genes = {}
    for x in new_sub_refbed:  
        #print(x)
        rg = [y - refshift for y in x[0:2]]
        h0,h1 = h_rbox
        ele = x[2]
        col = ele_col(ele)
        add_ele_box (ax, rg, h0, h1, fc=col, ec=col)
        gene_k = x[3] + ' ' + x[4]
        if gene_k not in ref_genes : ref_genes[gene_k] = []
        ref_genes[gene_k].append(rg[0])
        ref_genes[gene_k].append(rg[1])
    for gene_k in ref_genes:
        pos_arr = ref_genes[gene_k]
        rg = [min(pos_arr), max(pos_arr)]
        x = sum(rg)/2
        h0,h1 = h_rtext
        y = (h0 + h1)/2
        text = gene_k[3:].split()[0]
        size = 8
        if gene[3] == 'V' : 
            ha = 'left'
            x = rg[0]
        elif gene[3] == 'J' : 
            ha = 'right'
            x = rg[1]
        else : 
            ha = 'center'

        if check_overlap(rg, ref_intervals) :
            add_text(ax, text, x, y, size, ha=ha)
        h0,h1 = h_rbox
        add_mid_line(ax, rg, h0, h1, lw=2)

    ax.set_xticks([])
    ax.set_yticks([])


# load parameters



# load files
refbed, ref_ctg_list = load_bed(refbedfile)
tarbed, tar_ctg_list = load_bed(tarbedfile)
tarmeta = load_meta(tarmetafile)
refmeta = load_meta(refmetafile)
paf, keys = load_paf(paffile)
plt_nrow,plt_ncol = layout
fig_width,fig_height = figsize


if order_tag :
    tarmeta_order = load_meta_order(tarmetafile)
    new_keys = []
    for x in tarmeta_order :
        for i in range(len(keys)) :
            if x == keys[i].split()[0] :
                new_keys.append(keys[i])
    for x in keys :
        if x not  in new_keys:
            new_keys.append(x)
    keys = new_keys

print(keys)
#exit(0)

block_len = plt_nrow * plt_ncol
block_n = len(keys) // block_len + 1
key_set = []
for i in range(block_n) : key_set.append([])
for i in range(len(keys)) :
    idx = i // block_len
    key_set[idx].append(keys[i])

pdf_pages = PdfPages(outpref + '.pdf')
for ks in key_set :
    if not ks : break
    fig, axs = plt.subplots(plt_nrow, plt_ncol, 
                            figsize=(fig_width,fig_height),
                            squeeze=False)
    for ax in axs.ravel(): 
        ax.set_xticks([])
        ax.set_yticks([])

    plt_i = 0
    for k in ks :
        print('processing\t'  + k, file=sys.stderr)
        ctg,ref = k.split(' ')
        ctg0 = ctg.split(':')[0]
        tar_info = ''
        ref_info = ''
        if ctg0 in tarmeta : 
            tar_info = ';'.join(tarmeta[ctg0])
        elif ctg in tarmeta :
            tar_info = ';'.join(tarmeta[ctg])
        if ref in refmeta : ref_info = ';'.join(refmeta[ref])
        if ctg in tarbed :
            sub_tarbed = tarbed[ctg]
        else :
            sub_tarbed = []
        sub_refbed = refbed[ref]
        i = plt_i // plt_ncol
        j = plt_i % plt_ncol
        plot_one_pair(paf[k], axs[i][j], ctg, ref, sub_tarbed, sub_refbed, tar_info, ref_info, gapratio = gapratio)
        plt_i += 1
    plt.subplots_adjust(wspace=0, hspace=0)
    pdf_pages.savefig(bbox_inches='tight')
    plt.close()
pdf_pages.close()
