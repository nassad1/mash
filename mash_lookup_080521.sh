#!/bin/bash

#08/05/21: added background subsampling

FOREGROUND_PATH_MULTI=$1
BACKGROUND_PATH_MULTI=$2
METADATA_PATH=$3
#EXCLUDE_PATH=$4
MASH_PATH=$4
SEQUENCE_CLEANER_PATH=$5

module load gcc

#####BACKGROUND PROCESSING#####

local_path="$(dirname "${BACKGROUND_PATH_MULTI}")"
cd "$local_path"

mkdir background
local_path_background="${local_path}/background"
cd "$local_path_background"
cp $BACKGROUND_PATH_MULTI $local_path_background
BACKGROUND_PATH_MULTI_CP="${local_path_background}/$(basename "$BACKGROUND_PATH_MULTI")"

#remove sequences <27k length or >3,000N
ml python ####&> /dev/null
python $SEQUENCE_CLEANER_PATH $(basename "$BACKGROUND_PATH_MULTI_CP") 27000 3000 ####&> /dev/null

BACKGROUND_PATH_MULTI_CP_FILTERED="${local_path_background}/filtered_$(basename "$BACKGROUND_PATH_MULTI")"








###remove identical background sequences with same metadata
#echo "Checking for duplicate sequences."
#convert fasta to tsv
awk -v RS=">" -v ORS="\n" -v OFS="" '{$1="#"$1"\t"}1' $BACKGROUND_PATH_MULTI_CP_FILTERED|tail -n+2 > sorted_combined.fasta
#append metadata to sequence(if country is USA: "sequence"_"date"_"country"_"division"; if country is not USA: "sequence"_"date"_"country")
{
read
while IFS=$'\t' read -r p;
do
COUNTRY=$(echo "$p" | awk -F"\t" '{print $7}')
COUNTRY="$(echo $COUNTRY | sed 's/ /\__/g')"
if [ $COUNTRY = "USA" ]
then
   echo "$p" | awk -F"\t" '{print $5 "_" $7 "_" $8}' >> temp_metadata.tsv
else
   echo "$p" | awk -F"\t" '{print $5 "_" $7}' >> temp_metadata.tsv
fi
done
} < $METADATA_PATH

paste -d'\t' sorted_combined.fasta temp_metadata.tsv > temp_combined_seq_metadata.tsv
awk -F"\t" '{print $1}' temp_combined_seq_metadata.tsv > temp_samplename.tsv
awk -F"\t" '{print $2 "_" $3}' temp_combined_seq_metadata.tsv > temp_seq_metadata.tsv
paste -d'\t' temp_samplename.tsv temp_seq_metadata.tsv > combined_seq_metadata.tsv
sort -u -k2,2 combined_seq_metadata.tsv > unique_seq_metadata.tsv
NUM_UNIQ_HITS=$(wc -l  < unique_seq_metadata.tsv)
#echo "${NUM_UNIQ_HITS} unique sequences found."




##was this just separating individual background lists and appending foreground sequence?

#isolate unique sample name vector in original GISAID format
awk -F"\t" '{print $1}' unique_seq_metadata.tsv | cut -c 2- | sed 's/__/\//g' > unique_samplename.tsv
#filter multifasta for unique output
sed -i -- 's/__/\//g' $BACKGROUND_PATH_MULTI_CP_FILTERED
temp=${FOREGROUND_PATH##*/}
temp=${temp%.*}
fasta_name="closest_unique_${temp}.fasta"
metadata_name="closest_unique_metadata_${temp}.tsv"
##### include foreground######
cp $FOREGROUND_PATH $fasta_name
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' unique_samplename.tsv combined.fasta >> $fasta_name
#filter metadata for unique output
head -1 combined_metadata.tsv > $metadata_name
grep -w -F -f unique_samplename.tsv combined_metadata.tsv >> $metadata_name







#replace global background samplename including / with __ so can create individual files with samplename (note use two _ to make unique)
sed -i -- 's/\//__/g' $BACKGROUND_PATH_MULTI_CP_FILTERED

#split multifasta into individual files
pip install pyfasta ####&> /dev/null
pyfasta split --header "%(seqid)s.fasta" $BACKGROUND_PATH_MULTI_CP_FILTERED

rm $BACKGROUND_PATH_MULTI_CP
rm $BACKGROUND_PATH_MULTI_CP_FILTERED

#remove foreground sequences submitted to GISAID
#rm *-PV* ####&> /dev/null
#rm *USA__SC-20* ####&> /dev/null

#sed -i -- 's/\//__/g' $EXCLUDE_PATH
#while read p; do
#  if [[ $p != *"#"* ]] && [[ $p != "" ]];then
#    rm -f *"$p"* ####&> /dev/null
#  fi  
#done < $EXCLUDE_PATH

##break global into multiple folders (prevent system error 2*10^6 max)
dir_size=10000
dir_name="global"
n=$((`find . -maxdepth 1 -type f | wc -l`/$dir_size+1))

echo $n############

for i in `seq 1 $n`;
do
    mkdir -p "$dir_name$i";

echo "$dir_name$i"############

    find . -maxdepth 1 -type f | head -n $dir_size | xargs -i mv "{}" "${dir_name}${i}"
done

cd $local_path_background

num_fold=$(ls | wc -l)
for ((i=1;i<=$num_fold;i++));
do
 $MASH_PATH sketch -p 64 -o global_sketch_${i} "${local_path_background}"/"${dir_name}""${i}"/*.fasta ############&> /dev/null

 echo "${dir_name}""${i}"############

 ##paste iteratively and remove previous sketch.
 #if [ ${i} == 1 ]
 #then
 # $MASH_PATH paste global_sketch global_sketch_*.msh ############&> /dev/null
 # rm global_sketch_*.msh
 #else
 # $MASH_PATH paste global_sketch global_sketch global_sketch_*.msh ############&> /dev/null
 # rm global_sketch_*.msh
 #fi

done

$MASH_PATH paste global_sketch global_sketch_*.msh ############&> /dev/null
#rm global_sketch_*.msh

BACKGROUND_SKETCH_PATH="${local_path_background}/global_sketch.msh"

#####FOREGROUND#####

local_path="$(dirname "${FOREGROUND_PATH_MULTI}")"
cd "$local_path"

mkdir foreground
cd foreground
local_path_foreground="$(pwd)"
cp $FOREGROUND_PATH_MULTI $local_path_foreground
FOREGROUND_PATH_MULTI_CP="$(pwd)/$(basename "$FOREGROUND_PATH_MULTI")"

#replace foreground samplename including / with __ so can create individual files with samplename (note use two _ to make unique)
sed -i -- 's/\//__/g' $FOREGROUND_PATH_MULTI_CP

#split foreground multifasta into individual files
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < $FOREGROUND_PATH_MULTI_CP

rm $FOREGROUND_PATH_MULTI_CP


##APPROACH 1
#FOR EACH FILE IN local_path_foreground
for filename in *.fa*
do 
  FOREGROUND_PATH="$filename"

#sketch foreground
$MASH_PATH sketch -p 64 -o input_fasta $FOREGROUND_PATH ####&> /dev/null

#calculate distance matrix and sort

$MASH_PATH dist -p 64 $BACKGROUND_SKETCH_PATH input_fasta.msh > mtx_global_input_compiled.tab
#####SEGMENTATION FAULT
#maybe have to specify k or whatever? they look consistent: 1,000 hashes, 29903 length genome, 21 k-mer length
#issue is with paste of background sketch
#loop through each instead of pasting into one large?


sort -k 3,3 mtx_global_input_compiled.tab > mtx_sorted.tab

#find lowest distance value given threshold .0003
MIN_DIST=$(awk 'NR==1 && $3 < .0003 {print $3}' mtx_sorted.tab)

if [ -n "$MIN_DIST" ];then
#extract path to all minimum distance background sequences
awk -v min="$MIN_DIST" '$3==min{print $1}' mtx_sorted.tab > min_background_path.tab
NUM_HITS=$(wc -l  < min_background_path.tab)

#echo "${NUM_HITS} sequences found with minimum hash distance of ${MIN_DIST}."

#print foreground and non-identical lowest distance background to multi-fasta
#echo "Printing multi-fasta."
cp $FOREGROUND_PATH combined.fasta
temp_path=$(awk '{print $1}' min_background_path.tab)
cat $temp_path >> combined.fasta

#print list-view with metadata
#echo "Printing metadata."
#add header to label columns of combined_metadata.tsv
head -1 $METADATA_PATH > combined_metadata.tsv

#use minimum background sequence path to find corresponding metadata (convert filename back to sample name via __ to /)
while IFS=$'\t' read -r -a myArray
do
 temp=${myArray[0]}
 temp=${temp##*/}
 temp=${temp%.*}
 temp="$(echo $temp | sed 's/__/\//g')"
 awk -v sample="$temp" '$1==sample{print $0}' $METADATA_PATH >> combined_metadata.tsv
done < min_background_path.tab

#remove identical sequences with same metadata
#echo "Checking for duplicate sequences."
#convert fasta to tsv
awk -v RS=">" -v ORS="\n" -v OFS="" '{$1="#"$1"\t"}1' combined.fasta|tail -n+2 > sorted_combined.fasta
#remove foreground row
echo "$(tail -n +2 sorted_combined.fasta)" > sorted_combined.fasta
#append metadata to sequence(if country is USA: "sequence"_"date"_"country"_"division"; if country is not USA: "sequence"_"date"_"country")
{
read
while IFS=$'\t' read -r p;
do
COUNTRY=$(echo "$p" | awk -F"\t" '{print $7}')
COUNTRY="$(echo $COUNTRY | sed 's/ /\__/g')"
if [ $COUNTRY = "USA" ]
then
   echo "$p" | awk -F"\t" '{print $5 "_" $7 "_" $8}' >> temp_metadata.tsv
else
   echo "$p" | awk -F"\t" '{print $5 "_" $7}' >> temp_metadata.tsv
fi
done
} < combined_metadata.tsv

paste -d'\t' sorted_combined.fasta temp_metadata.tsv > temp_combined_seq_metadata.tsv
awk -F"\t" '{print $1}' temp_combined_seq_metadata.tsv > temp_samplename.tsv
awk -F"\t" '{print $2 "_" $3}' temp_combined_seq_metadata.tsv > temp_seq_metadata.tsv
paste -d'\t' temp_samplename.tsv temp_seq_metadata.tsv > combined_seq_metadata.tsv
sort -u -k2,2 combined_seq_metadata.tsv > unique_seq_metadata.tsv
NUM_UNIQ_HITS=$(wc -l  < unique_seq_metadata.tsv)
#echo "${NUM_UNIQ_HITS} unique sequences found."

#isolate unique sample name vector in original GISAID format
awk -F"\t" '{print $1}' unique_seq_metadata.tsv | cut -c 2- | sed 's/__/\//g' > unique_samplename.tsv
#filter multifasta for unique output
sed -i -- 's/__/\//g' combined.fasta
temp=${FOREGROUND_PATH##*/}
temp=${temp%.*}
fasta_name="closest_unique_${temp}.fasta"
metadata_name="closest_unique_metadata_${temp}.tsv"
##### include foreground######
cp $FOREGROUND_PATH $fasta_name
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' unique_samplename.tsv combined.fasta >> $fasta_name
#filter metadata for unique output
head -1 combined_metadata.tsv > $metadata_name
grep -w -F -f unique_samplename.tsv combined_metadata.tsv >> $metadata_name

#make output for IDs
ID_temp="ID_temp_${temp}.tsv"
ID_min_dist="ID_min_dist_${temp}.tsv"
ID_output="ID_output_${temp}.tsv"
foreground="$(echo $temp | sed 's/__/\//g')"
yes ${foreground} | head -n $NUM_UNIQ_HITS >> $ID_temp
yes ${MIN_DIST} | head -n $NUM_UNIQ_HITS >> $ID_min_dist
paste $ID_temp unique_samplename.tsv $ID_min_dist >> $ID_output

fi

rm input_fasta.msh ####&> /dev/null
rm mtx_global_input_compiled.tab ####&> /dev/null
rm mtx_sorted.tab ####&> /dev/null
rm min_background_path.tab ####&> /dev/null
rm combined.fasta ####&> /dev/null
rm combined_metadata.tsv ####&> /dev/null
rm sorted_combined.fasta ####&> /dev/null
rm temp_metadata.tsv ####&> /dev/null
rm temp_combined_seq_metadata.tsv ####&> /dev/null
rm temp_samplename.tsv ####&> /dev/null
rm temp_seq_metadata.tsv ####&> /dev/null
rm combined_seq_metadata.tsv ####&> /dev/null
rm unique_seq_metadata.tsv ####&> /dev/null
rm unique_samplename.tsv ####&> /dev/null
#rm $fasta_name ####&> /dev/null
#rm $metadata_name ####&> /dev/null
rm $ID_temp ####&> /dev/null
rm $ID_min_dist ####&> /dev/null

done


cat ID_output_*.tsv > master_output.tsv

cat ID_output_*.tsv
