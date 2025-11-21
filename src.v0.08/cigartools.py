import sys, re

def parsing_cs(cs) :
    matches = re.findall(r'([:\+\-\*]{1})([\da-z]+)', cs)
    out = [[m[0], m[1]] for m in matches]
    return(out)

def rev_complement(seq) :
    rc_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
              'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    rev_seq = seq[::-1]
    return(''.join([rc_map[x] for x in rev_seq]))



def flip_cs(cs) :
    # flip cigar string from query-to-ref to ref-to-query
    mat = [x for x in reversed(parsing_cigar(cgr))]
    out = []
    for i in range(len(mat)):
        if mat[i][0] == '+' : 
            mat[i][0] = '-'
            mat[i][1] = rev_complement(mat[i][1])
        elif mat[i][0] == '-' : 
            mat[i][0] = '+'
            mat[i][1] = rev_complement(mat[i][1])
        out.append(mat[i][0]+ mat[i][1])
    return(''.join(out))

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

def cal_match_score_cs(cs) :
    mat = parsing_cs(cs)
    score = 0
    for x in mat :
        tag = x[0]
        info = x[1]
        if tag == ':' : score += int(info)
        elif tag == '-' : score -= len(info)
        elif tag == '+' : score -= len(info)
        elif tag == '*' : score -= 1
        else :
            exit("unknow cs tag: " + tag)
    return(score)

def parsing_cigar(cgr) :
    matches = re.findall(r'(\d+)([A-Z=]{1})', cgr)
    out = [[m[0], m[1]] for m in matches]
    return(out)

def flip_cigar(cgr, strand) :
    # flip cigar string from query-to-ref to ref-to-query
    if strand == '-' :
        mat = [x for x in reversed(parsing_cigar(cgr))]
    else :
        mat = parsing_cigar(cgr)
    out = []
    for i in range(len(mat)):
        if mat[i][1] == 'D' : mat[i][1] = 'I'
        elif mat[i][1] == 'I' : mat[i][1] = 'D'
        out.append(mat[i][0]+mat[i][1])
    return(''.join(out))

def extract_cigar(fro, to, cgr) :
    # 0-based coordinate
    all_tab = parsing_cigar(cgr)
    new_tab = []
    cur_pos = 0
    cgr_l = 0
    for x in all_tab :
        l = int(x[0])
        tag = x[1]
        if tag in ['=','X','D'] : cgr_l += l
    for x in all_tab :
        l = int(x[0])
        tag = x[1]
        if tag in ['=','X','D'] : inc = l
        elif tag == 'I' : inc = 0
        else :
            exit("unknow cigar tag: " + tag)
        new_pos = cur_pos + inc
        # check orverlaping between [fro, to) and [cur_pos, new_pos)
        if cur_pos >= to : break
        elif new_pos <= fro : 
            pass
        else :
            new_l = min(to, new_pos) - max(fro, cur_pos)
            new_tab.append(str(new_l) +tag)
        #print("->", fro, to, cur_pos, new_pos, x, new_tab)
        cur_pos = new_pos
    return([''.join(new_tab), cgr_l])

def cal_match_score(fro, to, cgr) :
    cgr0, cgr_l = extract_cigar(fro, to, cgr)
    mat = parsing_cigar(cgr0)
    score = 0
    for x in mat :
        tag = x[1]
        l = int(x[0])
        if tag == '=' : score += l
        elif tag == 'X' : score -= 5
        elif tag == 'D' : score -= 10
        elif tag == 'I' : score -= 10
        else :
            exit("unknow cigar tag: " + tag)
    return(score)



if __name__ == '__main__':
    print("test CS string module:")
    cs = ':10+agg:5+a:2-t*ga:3-ccc*ca:1+a:2-ct:3+tatttt:1*cg:9'
    print("Input: " + cs)
    print("test parsing:")
    for x in parsing_cs(cs) :
        print(x)
    print("test extracting, 0-based coordinate")
    fro = 12; to = 20
    print(fro, to)
    subcs = extract_cs(fro, to, cs)
    print(subcs)
    print("match score of the extracted region:")
    print(cal_match_score_cs(subcs))



