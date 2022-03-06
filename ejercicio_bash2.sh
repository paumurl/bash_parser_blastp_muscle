#!/bin/bash
#Ejercicio Evaluable bash v1.0
#Paula Robles López, 2021

# This script allows to performs a blastp analysis using a urls_file to download the subject proteomes
# this file have to contain two columns, first column contains the species (identifier) for that proteome
# and the url. The script will generate a output project folder which contains two subfolders, data (here
# downloaded proteomes are stored and input fasta files) and results (here blast result, in addition a folder is created
# for every query protein to store blast result, aligment and trees for each protein)



#Function to print help in terminal
ayuda () {
   echo 'Ejercicio Evaluable bash v1.0'
   echo 'Paula Robles López, 2021'
   echo -e '\nusage: ejercicio.sh <query_sequences.fa> <ncbi_urls_file> <output_folder> <blast-identity> <blast-coverage>\n'
   echo -e 'query_sequences.fa : a fasta/multifasta file containing protein query sequences for blast'
   echo -e 'ncbi_urls_file     : a text plain file containing species name and url to download fasta protein file'
   echo -e 'output_folder      : folder in which data and results will be stored'
   echo -e 'blast-identity     : sequence identity cut off value 0-100'
   echo -e 'blast-coverage     : sequence coverage cut off value 0-100\n'
}


#Arguments assignation
query=$1
ncbi_urls_file=$2
project_name=$3
iden=$4 #70
cov=$5 #40


#Controls help message
echo -e "\n"
if (( $#==1 )); then #si solo hay un argumento
	case "$1" in
		"-h") ayuda; exit 0;; # y coincide que es petición de ayuda llamamos a la función de ayuda
		"-help") ayuda; exit 0;;
		*) echo -e "error:too few arguments"; ayuda; exit 0;; #si no tenemos suficientes argumentos
	esac
fi

#Control of arguments number 
if (( $#<5 )); then #mensaje de error si no tenemos suficientes argumentos
	echo -e "error:too few arguments"
	ayuda
	exit 0
fi




#Create project directories and output_files
if [[ -e $3 ]]; then #si ya existe un directorio con el nombre que nos da el usuario
	echo -e "Warning! That output folder already exists, please use another name\n"
	exit 0 #nos da error y salimos
else #si no pues creamos el proyecto y data y results dentro del nuevo directorio proyecto
	mkdir $3
	mkdir $3/data
	mkdir $3/results
fi


#parsing $ncbi_urls_file
echo -e "$0 is Running...\n"
echo "Fasta files will be stored in /$project_name/data/ "

cat $ncbi_urls_file | while read line #leemos el archivo con las urls del ncbi para descargar los .gz
do
	linea=$(echo $line | tr " " "\t" | cut -f2) #nos quedamos con el link de descarga
	echo "Downloading fasta files from NCBI..."
	wget -P $3/data $linea 1>>$3/log 2>>$3/log #descargamos cada .gz y guardamos en data
done

for archivo in $(ls $3/data/*.gz);do #descomprimimos todos los archivos que hemos descargado
	gzip -d $archivo
done

cat $3/data/*.faa >$3/data/proteome.fa #archivo fasta con todas las proteínas de todos los .faa que hemos descomprimido


cp $query $3/data #copiamos el archivo query a data


#Generate a big file that contains for every species all the fasta headers, necesary for species assgination in blast hits
echo "Generating species_fasta_ID.tsv file... "

for fasta in $(ls $3/data/*.faa);do #para cada archivo .faa que hemos descomprimido
	file=$(basename $fasta) #nos quedamos con el nombre del archivo solo y no toda su ruta
        nombre=$(grep "$file" $2 | tr " " "\t" | cut -f1) #nos quedamos con el nombre del organismo de cada .faa
	echo -e "$nombre\t$fasta" >> $3/data/Identifiers.tsv #imprimimos el organismo y la ruta de su proteoma
        grep "^>.*\s" $3/data/$file | tr " " "\t" | cut -f1 > $3/data/archivo_seq  # grepeamos los headers del .faa para ese organismo
        cat $3/data/archivo_seq | while read line; do #del archivo donde guardamos los headers leemos cada línea
                echo -e "$nombre\t$line" >> $3/data/species_fasta_ID.tsv # e imprimimos las correspondencias de los headers con su bicho
        done
done

rm $3/data/archivo_seq #eliminamos el archivo temporal que he creado para facilitar el parseo


#### BLAST SECTION ####

#Blast analysis and result filtering
echo "Running Blast analysis and result filtering... "

#hacemos blast con nuestro query y multifasta que hemos creado antes
(blastp -query $3/data/$query -subject $3/data/proteome.fa -outfmt "6 qseqid sseqid pident qcovs sseq" 2>>$3/log) > $3/results/blast_result

#con el archivo anterior resultado del blasta hacemos filtro para la identidad y coverage especificados
cat $3/results/blast_result | awk -v iden="$iden" -v cov="$cov" '{if ($3>=iden && $4>=cov) {print $0}}' > $3/results/blast_result_filtered 

#de este nuevo archivo leemos y nos quedamos con cada columna, podría haber cogido varias columnas a la vez pero he preferido hacerlo así
cat $3/results/blast_result_filtered | while read line; do #para cada linea del archivo filtrado
	id1=$(echo $line | tr " " "\t" | cut -f1) #nos quedamos con la susodicha columna
	id2=$(echo $line | tr " " "\t" | cut -f2)
	pident=$(echo $line | tr " " "\t" | cut -f3)
	qcovs=$(echo $line | tr " " "\t" | cut -f4)
	sseq=$(echo $line | tr " " "\t" | cut -f5)
#nos quedamos con el nombre del organismo que grepee en el tsv en el que tenemos cada header con su correspondiente organismo
	organismo=$(grep "$id2" $3/data/species_fasta_ID.tsv | tr " " "\t" | cut -f1)
	echo -e "$id1\t$organismo\t$id2\t$pident\t$qcovs\t$sseq">>$3/results/blast_result_final #nuevo archivo pero + el nombre del bicho

done

grep "^>.*" query.fa | sed 's/>//g' > $3/results/proteinasquery #nos quedamos con el número de secuencias query que tenemos

cat $3/results/proteinasquery | while read line; do #para cada secuencia query
	mkdir $3/results/$line #creamos una carpeta
	grep $line $3/results/blast_result_final > $3/results/$line/$line.1.fa #grepeamos la linea entera que sea hit para  la seq query
#formateo y seleccion de datos para obtener un multifasta con el organismo_fasta_header como nuevo header
	awk '{print ">"$2"_"$3"\n"$6}' $3/results/$line/$line.1.fa > $3/results/$line/$line.fa
	rm $3/results/$line/$line.1.fa #eliminamos un archivo temporal que he usado para que no hubiese errores en el volcado
done

### MUSCLE SECTION ####

#Make MUSCLE Phylogenetic trees for each query protein
echo "Making MUSCLE alignment and phylogenetic trees..."

cat $3/results/proteinasquery | while read line; do #para cada secuencia query
	muscle -in $3/results/$line/$line.fa -out $3/results/$line/$line.aln 1>>$3/log 2>>$3/log #hacemos muscle
	muscle -maketree -in $3/results/$line/$line.aln -out $3/results/$line/$line.nw -cluster neighborjoining 1>>$3/log 2>>$3/log
	#creamos el arbol con el resultado del alineamiento del muscle

done

### TERMINAL PRINTING ####

# Final stats to show in terminal
echo -e "\n\n"
echo -e "*** Done!! ***\n"
echo -e "Results are available at $3/results !"
echo -e "hits were found for query proteins:"
suma=0
#para cada secuencia query
cat $3/results/proteinasquery | \
{ while read line; do
        num=$(grep -c ">" $3/results/$line/$line.fa) #grepeamos cuantos hits hay para cada secuencia query
        echo -e -n "$num $line " #imprimimos num hits  y el nombre de la secuencia query
        suma=$((suma+num)) #sumamos los hits totales
done
echo -e "\nTotal blast hits: $suma\n" #hits totales
}
