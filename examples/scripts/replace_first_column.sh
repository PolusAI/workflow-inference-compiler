# See https://stackoverflow.com/questions/7846476/replace-column-in-one-file-with-column-from-another-using-awk

cat $1 | grep -v "#" | grep -v "@" > tempone.txt
cat $2 | grep -v "#" | grep -v "@" > temptwo.txt
awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' temptwo.txt tempone.txt