#/usr/bin/perl -w
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(13);
use Time::HiRes 'time';
#use warnings;
=pod
GetOptions( 'banco=s' => \$banco,
'organismo=s' => \$organismo,
'help=s' => \$help,
'output=s' => \$out,
);
=cut
########### VARIAVEIS DO SISTEMA ############ (NAO MODIFIQUE)
$consultaID = time();
################ PARAMETROS #################
#--------GERAL--------#
$dirTEMP = "/temp/"; #Pasta que sera usada pelo sistema para armazenar arquivos temporarios
$patogeno = "/home/arthur/databases/toxoplasma/sequencias_toxo.fasta"; #Proteoma do patogeno (organismo1)
$hospedeiro = "/home/tiago/bancos_de_dados/human.dir/Homo_sapiens.GRCh37.67.pep.all.fa"; #Proteoma do hospedeiro (organismo2)
#--------PEIMAP---------#
$basePEIMAP = "/home/tiago/peimap_scop100_nr70_20090102_no_string.mpfa"; #Proteoma do banco PEIMAP
$dirArquivosPEIMAP = "/home/thiago/arquivos_PEIMAP/PEIMAP_partes"; #Diretorio dos arquivos de interacoes PEIMAP
$arquivoReliabilidadeMetodos = "peimap/reliabilidade_metodos.txt"; #Arquivo que armazena os valores de reliabidade dos metodos experimentais
$numeroPEIMAPpartes = 10; #default = 3951
#---------PFAM----------#
$basePFAM = "/home/arthur/arquivos_Pfam/"; #Proteoma do banco PFAM (diretorio das matrizes .HMM)
$dirArquivosPFAM = "/home/arthur/arquivos_Pfam/iPfam_files/interactions_iPfam_triado_sem_redundancias.txt"; #Arquivos de interacoes PFAM
#--------PSIMAP---------#
$basePSIMAP = "/home/arthur/arquivos_SCOP/documentos_raiz/sequencias_SCOP_fasta/SCOP_isoformas_oficial.fasta"; #Proteoma do banco PSIMAP
$dirArquivosPSIMAP = "/home/arthur/arquivos_PSIMAP/PSIMAP_partes"; #Diretorio dos arquivos de interacoes PSIMAP
$numeroPSIMAPpartes = 10; #default = 598
###############################################

sub main(){
   print "Iniciando PIPELINE - Consulta n$consultaID\n";
   #$pm->start and next;
      runPEIMAP();
   #$pm->finish;
   #$pm->start and next;
      #runPFAM();
   #$pm->finish;
   #$pm->start and next;
      #runPSIMAP();
   #$pm->finish;
	
   
   
   #scoreFinal();	

}
main();

sub runPEIMAP(){
print "\tConsulta $consultaID:PEIMAP\n";
	##-----1º step
	#Run Blastp alignment between the two input proteins
   print "blast\n";
   print "blastall -p blastp -d $basePEIMAP -i $patogeno -F T -m8 -o "."PEIMAP$consultaID"."_pat1\n";
   print "----------blast\n";
	system("blastall -p blastp -d $basePEIMAP -i $patogeno -F T -m8 -o "."PEIMAP$consultaID"."_pat1");
   print "blastall -p blastp -d $basePEIMAP -i $hospedeiro -F T -m8 -o "."PEIMAP$consultaID"."_hos1\n";
   system("blastall -p blastp -d $basePEIMAP -i $hospedeiro -F T -m8 -o "."PEIMAP$consultaID"."_hos1");
   
   #blastall -p blastp -d /home/tiago/peimap_scop100_nr70_20090102_no_string.mpfa -a 10 -i /home/arthur/databases/toxoplasma/sequencias_toxo.fasta -F T -m8 -o out_blastp.txt
   #input blast organismo: /home/arthur/databases/toxoplasma/sequencias_toxo.fasta
   
	##--------2º step
	#Filter 30% < x < 70%
	system("PEIMAP$consultaID"."_pat1 | awk '$3 > 30 && $4 > 70{print $1\" \"$2}'  > PEIMAP$consultaID"."_pat2");
   system("PEIMAP$consultaID"."hos1 | awk '$3 > 30 && $4 > 70{print $1\" \"$2}'  > PEIMAP$consultaID"."_hos2");
   
	##-------3º step
	#Remove repetition with uniq command
	system("uniq PEIMAP$consultaID"."_pat2 > PEIMAP$consultaID"."_pat3");
   system("uniq PEIMAP$consultaID"."_hos2 > PEIMAP$consultaID"."_hos3");
   
	##-------4º step
	#Constructing the interaction network
	for($i=0;$i<=$numeroPEIMAPpartes;$i++){
      $pm->start and next;
      system("python peimap/PEIMAPpyrede.py $dirArquivosPEIMAP"."/peimap_$i $patogeno $hospedeiro $arquivoReliabilidadeMetodos $dirTEMP"."PEIMAP$consultaID"."_tmp$i");
      $pm->finish;
   }
	##Join all parts of the result into one file
	open(out,">PEIMAP$consultaID"."_pat+hos4");
	for($i=0;$i<=$numeroPEIMAPpartes;$i++){
		open(esp,$dirTEMP."PEIMAP$consultaID"."_tmp$i");
		@aux = <esp>;
		close(esp);
		foreach (@aux){
			print out"$_";
		}
	}
   close(out);

	##5º step retira redundancia e Juntas as duas primeiras colunas por um *, e elimina 3 e 4.
	system("sort -u PEIMAP$consultaID"."_pat+hos4 | awk '{print $1\"*\"$2\"\t\"$5}' > PEIMAP$consultaID"."_pat+hos5");
   system("sort -u PEIMAP$consultaID"."_pat+hos4 | awk '{print $1\"*\"$2}' > PEIMAP$consultaID"."_pat+hos5n");
   
   ##6º step
	#Calculate the n-value and normalize the reliability score
   system("uniq -c PEIMAP$consultaID"."_pat+hos5n > PEIMAP$consultaID"."_pat+hosNvalue");
	system("perl peimap/somador_de_reliabilidade_de_varias_interacoes.pl PEIMAP$consultaID"."_pat+hos5 > PEIMAP$consultaID"."_pat+hos6");
   
   ##7º step
   #Criar script para encontrar maior valor e normalizar todos os valores listados e unir os arquivos TGE*ENS scorenormalizado n-value
	system("python peimap/normalizador_PEI.py PEIMAP$consultaID"."_pat+hos6 PEIMAP$consultaID"."_pat+hosNvalue > PEIMAP$consultaID"."_pat+hos7");
}

sub runPFAM(){
  print "\tConsulta $consultaID:PFAM\n";
   ##1º Step
   system("perl /home/arthur/arquivos_Pfam/PfamScan/pfam_scan.pl -e_seq 0.05 -fasta $patogeno -dir /home/arthur/arquivos_Pfam/ -outfile "."PFAM$consultaID"."_pat1");
   system("perl /home/arthur/arquivos_Pfam/PfamScan/pfam_scan.pl -e_seq 0.05 -fasta $hospedeiro -dir /home/arthur/arquivos_Pfam/ -outfile "."PFAM$consultaID"."_hos1");
   #Ex: perl pfam_scan.pl -e_seq 0.05 -fasta /home/arthur/databases/toxoplasma/sequencias_toxo.fasta -dir /home/arthur/arquivos_Pfam/ -outfile pfam_out1.txt
   

   #Quando sair o resultado do pfam_scan, deve-se remover o cabeçalho
   ##2º Step
   #Extract the columns number 1, 6, 12
   system("cat PFAM$consultaID"."_pat1 | awk '{print $1\"\t\"$6\"\t\"$12}' > PFAM$consultaID"."_pat2");
   system("cat PFAM$consultaID"."_hos1 | awk '{print $1\"\t\"$6\"\t\"$12}' > PFAM$consultaID"."_hos2");
   
   ##3º Step
   #OBS: Além disso, o código Pfam permanecerá apenas com os  dígitos antes do ponto, para que este se relacione com as entradas do iPfam – banco que armazena as interações entre códigos Pfam.
   system("cat PFAM$consultaID"."_pat2 | awk -F \".\" '{print $1\"\t\"$2\".\"$3}' | awk '{print $1\"\t\"$2\"\t\"$4}' > PFAM$consultaID"."_pat3");
   system("cat PFAM$consultaID"."_hos2 | awk -F \".\" '{print $1\"\t\"$2\".\"$3}' | awk '{print $1\"\t\"$2\"\t\"$4}' > PFAM$consultaID"."_hos3");
   
   ##4º Step
   #Eliminate duplicated lines com as inversoes
=pod
   cat pfam_out3 | awk -F "\t" '{if($2 < $1){print $2"\t"$1"\t"$3}else{print $0}}' | sort -u | uniq > pfam_out4
   system("cat pfam_out3 | awk -F "\t" '{if($2 < $1){print $2"\t"$1"\t"$3}else{print $0}}' | sort -u | uniq > pfam_out4"); #tem que testar
   #Execute uniq command to elimitate duplicate records
	system("cat interespecifico_pre-uniq2 | uniq -f 5 -f 6 -f 7 -f 8 > interespecifico_triado");
=cut   

	##5º Step
	##Constructing the interaction network
	#system("perl interatoma_interespecifico_Pfam.pl /home/arthur/arquivos_Pfam/iPfam_files/interactions_iPfam_triado_sem_redundancias.txt $org.Pfamout $hum.Pfamout pfam_out4.txt");
   system("python pfam/PFAMpyrede.py $basePFAM $patogeno $hospedeiro PFAM$consultaID"."_pat+hos4");
	
	##6º Step
	#Eliminate duplicated lines
	system("sort -u PFAM$consultaID"."_pat+hos4 > PFAM$consultaID"."_pat+hos5");
   
   ##7º
   #Seleciona maior score
   system("perl pfam/seleciona_maior_score_Pfam.pl PFAM$consultaID"."_pat+hos5 > PFAM$consultaID"."_pat+hos6");
   
   #8º
   #Separa as duas primeiras colunas TGE*ENS e calcula n-vale
   system("cat PFAM$consultaID"."_pat+hos6 | awk '{print $1}' | uniq -c > PFAM$consultaID"."_pat+hosNvalue");
   
   #9º step
   #Criar script para encontrar maior valor e normalizar todos os valores listados e unir os arquivos TGE*ENS scorenormalizado n-value
   system("python pfam/normalizador_PF.py PFAM$consultaID"."_pat+hos6 PFAM$consultaID"."_pat+hosNvalue > PFAM$consultaID"."_pat+hos7");
}

sub runPSIMAP(){
   print "\tConsulta $consultaID:PSIMAP\n";
	system("blastpgp -d basePSIMAP -i $patogeno -F T -s T -j 5 -m8 > PSIMAP$consultaID"."_pat1");
   system("blastpgp -d basePSIMAP -i $hospedeiro -F T -s T -j 5 -m8 > PSIMAP$consultaID"."_hos1");
   
   system("cat PSIMAP$consultaID"."_pat1 | awk '{if($11 <= 0.0001){print $1\"\t\"$2}}' | sort -u > PSIMAP$consultaID"."_pat2");
   system("cat PSIMAP$consultaID"."_hos1 | awk '{if($11 <= 0.0001){print $1\"\t\"$2}}' | sort -u > PSIMAP$consultaID"."_hos2");
   
   
	#Constructing the interaction network
	
   for (my $i = 0 ; $i <= $numeroPSIMAPpartes; $i++) {
	$pm->start and next;
      system("python psimap/PSIMAPpyrede.py $dirArquivosPSIMAP/scop_$i $patogeno $hospedeiro $dirTEMP"."PSIMAP$consultaID"."_tmp$i");
	$pm->finish;
   }
	##Join all parts of the result into one file
	open(out,">PSIMAP$consultaID"."_pat+hos3");
	for($i=1;$i<$numeroPSIMAPpartes;$i++){
		open(esp,"$dirTEMP"."PSIMAP$consultaID"."_tmp$i");
		@aux = <esp>;
		close(esp);
		foreach (@aux){
			print out"$_";
		}
	}
	close(out);
	
	#aparentemente a menor distancia jA e selecionada no interatoma Teste e comparar depois
	#system("perl selecionador_de_pares_com_menor_distancia_PSIMAP.pl psimap_out4.txt > psimap_out5.txt");
	
	##Obtem-se as duas primeiras colunas e conta quantas vezes uma interação ocorreu atraves do uniq -c
	system("cat PSIMAP$consultaID"."_pat+hos3 | awk '{split($0,a,"*");print a[1]"*"a[2]}' | uniq -c > PSIMAP$consultaID"."_pat+hosNvalue");
   
   ##Separa as duas colunas TGE*ENS e score
   system("cat PSIMAP$consultaID"."_pat+hos3 | awk '{split($0,a,\"*\");print a[1]\"*\"a[2]\"\t\"$2}' > PSIMAP$consultaID"."_pat+hos4");
   
	#Criar script para encontrar maior valor e normalizar todos os valores listados e unir os arquivos TGE*ENS scorenormalizado n-value
   system("python psimap/normalizador_PSI.py PSIMAP$consultaID"."_pat+hos4 PSIMAP$consultaID"."_pat+hosNvalue > PSIMAP$consultaID"."_pat+hos5");
}

sub scoreFinal(){

}