import sys

filename = sys.argv[1] #interacoes preditas
filenvalue = sys.argv[2] #arquivo n-value
fnv = open(filenvalue)
nvalue = fnv.readlines()
fnv.close()
Dict = {}
for e in nvalue:
   col = e.split('\t')
   Dict[col[1]] = col[0]


file = open(filename)
val = []
linhas = []
for line in file:
   col = line.split('\t')
   linhas.append(col)
   val.append(float(col[1]))
maximo = max(val)
file.close()
####Normalizar
for i in range(len(linhas)):
   #linhas[i][1] = str(float(linhas[i][1])/maximo)
   print linhas[i][0]+'\t'+str(float(linhas[i][1])/maximo)+'\t'+Dict[linhas[i][0]]
