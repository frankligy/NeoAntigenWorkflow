

cat ineo_testing_filter910_new.txt | awk '$2=="HLA-A*6801"&&length($1)==9 {print $0}' > tmp7_9.txt
cat tmp7_9.txt | awk '{printf ">%s\n%s\n",$1,$1}' > tmp7_9.fasta
cat tmp6_9.txt | cut -f 3

# HLA-A0201
# HLA-A0101
# HLA-A1101
# HLA-A3101
# HLA-B0801
# HLA-B5801
# HLA-A6801