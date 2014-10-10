import sys

baseFilename = sys.argv[1]
patogenFilename = sys.argv[2]
hostFilename = sys.argv[3]
reliabilidadeFilename = sys.argv[4]
outputFilename = sys.argv[5]
#print sys.argv

file = open(baseFilename)
base = file.readlines()
file.close()

file = open(patogenFilename)
org = file.readlines()
file.close()

file = open(reliabilidadeFilename)
realibMetodos = file.readlines()
file.close()

file = open(hostFilename)
host = file.readlines()
file.close()



file = open(outputFilename,'w') #arquivo de saida
baseDict = {}
for e in base: #transforma os pares em hash, col1 eh chave para col2 e vice-versa
   col = e.split('\t')
   #print 'adicionando:', col
   if col[0] in baseDict: #chave ja existente
      if col[1] in baseDict[col[0]]: 
         baseDict[col[0]][col[1]] += [col[2].strip()]
      else:
         baseDict[col[0]][col[1]] = [col[2].strip()]
   else: #nao existe, temos de criar a chave
      baseDict[col[0]] = {col[1]:[col[2].strip()]}
   ##Fazer as chaves inversas
   if col[1] in baseDict: #chave ja existente
      if col[0] in baseDict[col[1]]: 
         baseDict[col[1]][col[0]] += [col[2].strip()]
      else:
         baseDict[col[1]][col[0]] = [col[2].strip()]
   else: #nao existe, temos de criar a chave
      baseDict[col[1]] = {col[0]:[col[2].strip()]}

#for k in baseDict.keys() : print k,baseDict[k]
#print baseDict
reliabDict = {}
for e in realibMetodos:
   col = e.split('\t')
   reliabDict[col[0]] = col[1]

for p in host:
   col = p.split('\t')
   print col
   if len(col) != 2 : continue
   col[1] = col[1].strip()
   print col
   if col[1] in baseDict: #verifica se existe interacao neste arquivo do banco
      pares = baseDict[col[1]].keys()
      #Encontramos o par, agora devemos procurar no organismo
      for o in org:
         coluna = o.split('\t')
         if len(coluna) != 2 : continue
         coluna[1] = coluna[1].strip()
         if coluna[1] in pares:
            metodos = min(baseDict[col[1]][coluna[1]]) #caracteristica PEIMAP
            metodos = metodos.replace(";",",")
            metodos = metodos.split(',')
            reliabilidade = 0.0
            for m in metodos:
               #print m, reliabDict[m]
               reliabilidade += float(reliabDict[m])
            ##Somar metodos
            file.write(col[0]+'\t'+coluna[0]+'\t'+col[1]+'\t'+coluna[1]+'\t'+str(reliabilidade)+'\n')
            #print col[0]+'\t'+coluna[0]+'\t'+col[1]+'\t'+coluna[1]+'\t'+str(reliabilidade)+'\n'
   
file.close()