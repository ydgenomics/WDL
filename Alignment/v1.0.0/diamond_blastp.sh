### Date: 250828
### Image: Alignment
### Ref: [DIAMOND:快又准的蛋白序列比对软件](https://mp.weixin.qq.com/s/5UhthY9PHfN7zxZbJdZaJA)
fasta1=$1
fasta2=$2
method=$3
n_cpu=$4
name1=$(basename "$fasta1")
name2=$(basename "$fasta2")

source /opt/software/miniconda3/bin/activate
conda activate alignment

mkdir result
if [[ "$method" == "diamond" ]]; then
  echo "Running diamond branch"
  diamond makedb --in $fasta1 --db $name1
  diamond makedb --in $fasta2 --db $name2
  diamond blastp --db $name2 -q $fasta1 -o "./result/blastp_"$name1"_vs_"$name2".txt"
  diamond blastp --db $name1 -q $fasta2 -o "./result/blastp_"$name2"_vs_"$name1".txt"
else
  echo "Running non-diamond branch"
  makeblastdb -in $fasta1 -dbtype prot -out $name1
  makeblastdb -in $fasta2 -dbtype prot -out $name2
  blastp -query $fasta1 -db $name2 -out "./result/blastp_"$name1"_vs_"$name2".txt" -outfmt 6 -evalue 1e-5 -num_threads $n_cpu
  blastp -query $fasta2 -db $name1 -out "./result/blastp_"$name2"_vs_"$name1".txt" -outfmt 6 -evalue 1e-5 -num_threads $n_cpu
fi

# get reciprocal result
echo -e "Query_ID\tRefer_ID\tIdentity(%)\tAlignment_Length\tMismatches\tGap_Openings\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score" > header.tsv
n=0
for i in $(ls */*.txt)
do 
  cat header.tsv $i > "$n".txt
  awk -F '\t' '$3 >= 70' "$n".txt > "$n"_filter.txt
  awk '!seen[$1]++' "$n"_filter.txt > "$n"_unique.tsv
  let n++
done

awk 'NR==FNR{a[$2"_"$1]=$1}NR!=FNR{if(a[$1"_"$2])print $1"\t"a[$1"_"$2]}' 0_unique.tsv 1_unique.tsv > reciprocal_best.txt