#########
## Set paths
input_read_dir="/path/to/dir"
output_dir="/path/to/out_dir"
sample_list_file="/path/to/file"

## Set parameters
threads=20
adap_fwd="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
adap_rev="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"
nruns=30

## Trim reads with AdapterRemoval
mkdir -p $output_dir/trimmed_reads
while read sample
do
  AdapterRemoval --identify-adapters --file1 $input_read_dir/${sample}*/*1.fq.gz --file2 $input_read_dir/${sample}*/*1.fq.gz --adapter1 $adap_fwd --adapter2 $adap_rev ---threads $threads --basename $output_dir/trimmed_reads/$sample
done < $sample_list_file

## Assemble genome
get_organelle_config.py --add all
while read sample
do
  mkdir -p $output_dir/${sample}_output
  get_organelle_from_reads.py -1 $output_dir/trimmed_reads/${sample}.pair1.truncated -2 $output_dir/trimmed_reads/${sample}.pair2.truncated -o $output_dir/${sample}_output -k 21,45,65,85,105 -F embplant_pt -R $nruns -t $threads
done < $sample_list_file

## Scaffold and circular genomes were annotated externally using GeSeq from Chlorobox. Markers along with upstream and downstream sequences were extracted using the respective annotation files.
# For each marker, download the gff3 annotation file with fasta sequences attached
marker="rbcL" #rbcL, matK, ndhF, trnF
while read sample
do
  fastastart=$(cat $sample.gff3 | grep -n ">" | cut -d':' -f1 | head -n1)
  tail -n +"$fastastart" $sample.gff3 > $sample.fasta
  if grep -q $marker $sample.gff3
  then
    data=$(cat $sample.gff3 | grep -P "\tgene\t" | grep $marker | awk '{print $5 - $4 "\t" $0}' | sort -k1,1 -rn | head -n1 | sed 's/;/\t/' | cut -f2,5,6,10) # Get longest sequence matching the marker name
    start=$(echo "$a" | cut -f2) # Start coordinate
    end=$(echo "$a" | cut -f3) # End coordinate
    sequence=$(echo "$a" | cut -f1) # Sequence name
    max=$(grep $sequence $sample.gff3 | grep -P "\tsource" | cut -f5) # Maximum coordinate of scaffold
    if [ "$start" -lt 5001 ]
    then
      let startext=1
    else
      let startext=$start-5000
    fi
    if [ "$end" -gt "$max" ]
    then
      let endext=$max
    else
      let endext=$end+5000
    fi
    name=$(echo "$data" | cut -f4 | grep -o $marker | sed "s/^/${sample}_/")
    seqkit subseq --chr $sequence -r ${startext}:$endext $sample.fasta > ${marker}_${sample}_5000ups_5000down.fasta
  else
    echo "Marker $marker not found for sample $sample"
  fi
done < $sample_list_file
