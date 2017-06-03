#if [ "$#" -lt 2 ]; then
#        echo "Usage: sh testpool.sh 10 0.5"
#        exit 1
#fi

for fname in $(ls *.fasta)
do
	echo $fname" "$(cat $fname | grep -o . |sort |uniq -c|wc -l)
done
