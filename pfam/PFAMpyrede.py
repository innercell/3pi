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
      baseDict[col[0]] += [col[1].strip()]      
   else: #nao existe, temos de criar a chave
      baseDict[col[0]] = [col[1].strip()]
   ##Fazer as chaves inversas
   if col[1] in baseDict: #chave ja existente
      baseDict[col[1]] += [col[0].strip()]      
   else: #nao existe, temos de criar a chave
      baseDict[col[1]] = [col[0].strip()]

#for k in baseDict.keys() : print k
#print baseDict

for p in host:
   col = p.split('\t')
   if col[0] == '#' or col[0] == '.' or len(col) != 3 : continue
   #col[1] = col[1].strip()
   #print col
   if col[1] in baseDict: #verifica se existe interacao neste arquivo do banco
      pares = baseDict[col[1]]
      #Encontramos o par, agora devemos procurar no organismo
      for o in org:
         coluna = o.split('\t')
         if coluna[0] == '#' or coluna[0] == '.' or len(coluna) != 3 : continue
         #coluna[1] = coluna[1].strip()
         if coluna[1] in pares:
            semelhanca1 = col[2]
            semelhanca2 = coluna[2]
            ### arquivo sem soma
            #file.write(col[0]+'\t'+col[2].strip()+'\t'+coluna[0]+'\t'+coluna[2].strip()+'\t'+col[1]+'\t'+coluna[1]+'\n')
            ### arquivo com soma das similaridades
            #file.write(col[0]+'*'+coluna[0]+'*'+col[1]+'*'+coluna[1]+'\t'+col[2].strip()+'\t'+coluna[2].strip()+'\t'+str(float(col[2].strip())+float(coluna[2].strip()))+'\n')
            file.write(col[0]+'*'+coluna[0]+'*'+col[1]+'*'+coluna[1]+'\t'+str(float(col[2].strip())+float(coluna[2].strip()))+'\n')
            
            #print col[0]+'\t'+col[2]+'\t'+coluna[0]+'\t'+coluna[2]+'\n'
   
file.close()