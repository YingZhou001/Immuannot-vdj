import sys


da = {}


seq = ''

for line in sys.stdin:
    if not line: break
    if line[0] == '>' :
        if seq :
            if seq not in da : da[seq] = []
            da[seq].append(metaline)
            seq = ''
            
        metaline = line[:-1]
    else :
        seq += line[:-1]

if seq not in da : da[seq] = []
da[seq].append(metaline)


for seq in da:
    print(da[seq][0])
    print(seq)
