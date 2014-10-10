$dirArquivosPEIMAP="/arquivos_PEIMAP/";
$patogeno = "peimap_org_out3.txt";
$hospedeiro = "peimap_host_out3.tx";
$arquivoReliabilidadeMetodos ="reliabilidade_metodos_PEIMAP.txt";
$dirTEMP= "/peimap_out/";
$consultaID="777";
##-------4º step
	#Constructing the interaction network
	for($i=0;$i<=3951;$i++){
      #print "parte $i";
      system("python PEIMAPpyrede.py $dirArquivosPEIMAP"+"peimap_$i $patogeno $hospedeiro $arquivoReliabilidadeMetodos $dirTEMP"+"PEIMAP$consultaID"+"_tmp$i");
   }
	##Join all parts of the result into one file
	open(out,">PEIMAP$consultaID"+"_pat+hos4");
	for($i=0;$i<=3951;$i++){
		open(esp,$dirTEMP+"PEIMAP$consultaID"+"_tmp$i");
		@aux = <esp>;
		close(esp);
		foreach (@aux){
			print out"$_";
		}
	}
   close(out);