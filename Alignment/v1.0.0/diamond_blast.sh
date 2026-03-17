### Date: 260317
### Image: Alignment
### Ref: [DIAMOND:快又准的蛋白序列比对软件](https://mp.weixin.qq.com/s/5UhthY9PHfN7zxZbJdZaJA)
fasta1=$1
fasta2=$2
type=$3 # nucleotide or protein
method=$4 # diamond or blast
n_cpu=$5
name1=$(basename "$fasta1")
name2=$(basename "$fasta2")

# # 使用grep检查是否存在点（只检查非标题行）
# if grep -v "^>" "$fasta1" | grep -q "\."; then
#     echo "发现序列中存在点(.)，正在替换为星号(*)..."
#     # 统计替换数量
#     dots_count=$(grep -v "^>" "$fasta1" | grep -o "\." | wc -l)
#     # 执行替换
#     sed '/^>/! s/\./\*/g' "$fasta1" > "$fasta1"
#     echo "完成！共替换了 $dots_count 个点"
# else
#     echo "序列中不存在点(.)，无需替换"
# fi

# # 使用grep检查是否存在点（只检查非标题行）
# if grep -v "^>" "$fasta2" | grep -q "\."; then
#     echo "发现序列中存在点(.)，正在替换为星号(*)..."
#     # 统计替换数量
#     dots_count=$(grep -v "^>" "$fasta2" | grep -o "\." | wc -l)
#     # 执行替换
#     sed '/^>/! s/\./\*/g' "$fasta2" > "$fasta2"
#     echo "完成！共替换了 $dots_count 个点"
# else
#     echo "序列中不存在点(.)，无需替换"
# fi

# source /opt/software/miniconda3/bin/activate && conda activate alignment

mkdir result
if [[ "$type" == "nucleotide" ]]; then
  echo "Use blastn alignement nucleotide sequences..."
  /opt/software/miniconda3/envs/alignment/bin/makeblastdb -in $fasta1 -dbtype nucl -out $name1
  /opt/software/miniconda3/envs/alignment/bin/makeblastdb -in $fasta2 -dbtype nucl -out $name2
  /opt/software/miniconda3/envs/alignment/bin/blastn -query $fasta1 -db $name2 -out "./result/blastn_"$name1"_vs_"$name2".txt" -outfmt 6 -evalue 1e-10 -num_threads $n_cpu
  /opt/software/miniconda3/envs/alignment/bin/blastn -query $fasta2 -db $name1 -out "./result/blastn_"$name2"_vs_"$name1".txt" -outfmt 6 -evalue 1e-10 -num_threads $n_cpu
else
  if [[ "$method" == "diamond" ]]; then
    echo "Use diamond alignement protein sequences..."
    /opt/software/miniconda3/envs/alignment/bin/diamond makedb --in $fasta1 --db $name1
    /opt/software/miniconda3/envs/alignment/bin/diamond makedb --in $fasta2 --db $name2
    /opt/software/miniconda3/envs/alignment/bin/diamond blastp --db $name2 -q $fasta1 -o "./result/blastp_"$name1"_vs_"$name2".txt"
    /opt/software/miniconda3/envs/alignment/bin/diamond blastp --db $name1 -q $fasta2 -o "./result/blastp_"$name2"_vs_"$name1".txt"
  else
    echo "Use blastp alignement protein sequences..."
    /opt/software/miniconda3/envs/alignment/bin/makeblastdb -in $fasta1 -dbtype prot -out $name1
    /opt/software/miniconda3/envs/alignment/bin/makeblastdb -in $fasta2 -dbtype prot -out $name2
    /opt/software/miniconda3/envs/alignment/bin/blastp -query $fasta1 -db $name2 -out "./result/blastp_"$name1"_vs_"$name2".txt" -outfmt 6 -evalue 1e-10 -num_threads $n_cpu
    /opt/software/miniconda3/envs/alignment/bin/blastp -query $fasta2 -db $name1 -out "./result/blastp_"$name2"_vs_"$name1".txt" -outfmt 6 -evalue 1e-10 -num_threads $n_cpu
  fi
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