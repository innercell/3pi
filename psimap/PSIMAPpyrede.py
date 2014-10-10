import sys

baseFilename = sys.argv[1]
patogenFilename = sys.argv[2]
hostFilename = sys.argv[3]
outputFilename = sys.argv[4]
#print sys.argv

file = open(baseFilename)
base = file.readlines()
file.close()

file = open(patogenFilename)
org = file.readlines()
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

#for k in baseDict.keys() : print k
#print baseDict

for p in host:
   col = p.split('\t')
   if len(col) != 2 : continue
   col[1] = col[1].strip()
   if col[1] in baseDict: #verifica se existe interacao neste arquivo do banco
      pares = baseDict[col[1]].keys()
      #Encontramos o par, agora devemos procurar no organismo
      for o in org:
         coluna = o.split('\t')
         if len(coluna) != 2 : continue
         coluna[1] = coluna[1].strip()
         if coluna[1] in pares:
            distMinima = min(baseDict[col[1]][coluna[1]]) #caracteristica PSIMAP
            if distMinima <= 0 : continue
            file.write(col[0]+'*'+coluna[0]+'*'+col[1]+'*'+coluna[1]+'\t'+distMinima+'\n')
            #print col[0]+'\t'+coluna[0]+'\t'+distMinima
   
file.close()