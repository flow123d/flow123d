import sys

def safediv(x, y):
    if y == 0:
        return 0
    else:
        return x / y

with open(sys.argv[1]) as source:
    content = source.readlines()
for i in range(20):
    print(content[i].strip('\n'))
print('%s %s %s %s' % ('2mw/w + mr/r', '     DLmw/Dw', '     DLmr/Dr', content[20].strip('\n')))
print(content[21].strip('\n'))
out = []
for i in range(22, len(content) - 1):
    columns = content[i].split()
    ratio1 = safediv(float(columns[8].replace(',', '')), float(columns[6].replace(',', '')))
    ratio2 = safediv(float(columns[5].replace(',', '')), float(columns[3].replace(',', '')))
    ratio3 = 2 * ratio1 + ratio2
    out.append([ratio3, '%.10f %.10f %.10f %s' % (ratio3, ratio1, ratio2, content[i].strip('\n'))])
out.sort(key=lambda x: x[0], reverse=True)
for line in out:
    print(line[1])
