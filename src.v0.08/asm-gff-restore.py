import sys

for line in sys.stdin:
    if not line: break
    if line[0] == '#' : 
        print(line, end='')
        continue
    llst = line[:-1].split('\t')
    if ':' in llst[0] and '-' in llst[0] :
        tmp = llst[0].split(':')
        fro,to = tmp[-1].split('-')
        ctg = ':'.join(tmp[:-1])
        llst[0] = ctg
        llst[3] = str(int(llst[3]) + int(fro) -1)
        llst[4] = str(int(llst[4]) + int(fro) -1)
    print('\t'.join(llst))
