fasta1="/data/input/Files/zhangzijian/cross/pep/Osativa_323_v7.0.protein1.fa"
fasta2="/data/input/Files/zhangzijian/cross/pep/Sbicolor_454_v3.1.1.protein1.fa"
method="blastp"
sh /WDL/Alignment/v1.0.0/diamond_blastp.sh $fasta1 $fasta2 $method 2