#!/usr/bin/env bash



echo "********Spoligotyping Started*********"
forReads=`ls | grep _R1`
echo "Forward Reads to be typed: $forReads"

spoligobinary=""
spoligocounts=""
spoligooctal=""

#94 Spacer Sequences


SS1="TTAAAACCGTGTTGCACTGCAACCCGGAATTCTTGCAC"
SS2="TTAC"
SS3="CAT"
SS4="GAT"



for i in {1..3}; do
	curseq=SS$i
	echo "SS to find $curseq, ${!curseq}, $curvar"
	seqcount=`grep -c ${!curseq} $forReads`
	echo "search on $i"
	spoligocounts="$spoligocounts $seqcount"
	done
	
echo "Spoligo Counts = $spoligocounts"

spoligobinary=`echo $spoligocounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
#brucbinary=`echo $bruccounts | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`

echo "$spoligobinary"

spoligooctal=`echo "ibase=2;obase=8; $spoligobinary" | bc`
#echo "ibase=2;obase=8; 1010101" | bc

echo "Spoligo Octal = $spoligooctal"

echo "Spoligo Counts = $spoligocounts" >> spoligotype.txt
echo "Spoligo Binary = $spoligobinary" >> spoligotype.txt
	

	






