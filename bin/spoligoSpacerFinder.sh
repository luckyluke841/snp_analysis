#!/bin/sh

#################################################################################
#  Dependencies ---
#   Unix..ish setup :)
#################################################################################

echo "**********************START**********************"

# Starting working directory must be BWA-GATK folder with included /zips file containing 2 zipped fastq files.
echo "directory"
pwd
# Make fastq directory
mkdir ./../fastq
cp ./../zips/*.fastq.gz ./../fastq

# change working directory to /fastq
cd ./../fastq

echo "starting to unzip files"
# Unzip files
gunzip *.fastq.gz
echo "finished unzipping files"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'` #grab name minus the .vcf

echo "2 3 4 12 13 14 15 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 46 47 51 52 53 62 63 64 65" > NUM

##############

onemismatch () {
patt=(${!forward})
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" > searchpattern

patt=(${!reverse})
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" >> searchpattern
}

##############
#Spacers- Original No./Extended No.

#Spacer -/01
forward1="CAACCCGGAATTCTTGC"
reverse1="GCAAGAATTCCGGGTTG"

#Spacer01/02
forward2="TGATCCAGAGCCGGCGACCCTCTAT"
reverse2="ATAGAGGGTCGCCGGCTCTGGATCA"

#Spacer02/03
forward3="CAAAAGCTGTCGCCCAAGCATGAGG"
reverse3="CCTCATGCTTGGGCGACAGCTTTTG"

#Spacer03/04
forward4="CCGTGCTTCCAGTGATCGCCTTCTA"
reverse4="TAGAAGGCGATCACTGGAAGCACGG"

#Spacer -/5
forward5="CAGGCGTGGCTAGG"
reverse5="CCTAGCCACGCCTG"

#Spacer -/6
forward6="GTCGCCGTAAGTGCC"
reverse6="GGCACTTACGGCGAC"

#Spacer -/7
forward7="GTTGACCACGAATTTTCAGA"
reverse7="TCTGAAAATTCGTGGTCAAC"

#Spacer -/8
forward8="GCTGGCGCGCATCAT"
reverse8="ATGATGCGCGCCAGC"

#Spacer -/9
forward9="CCATATCGGGGACGG"
reverse9="CCGTCCCCGATATGG"

#Spacer -/10
forward10="GCGTCGTGCCATCAG"
reverse10="CTGATGGCACGACGC"

#Spacer -/11
forward11="CCGTGCACATGCCGT"
reverse11="ACGGCATGTGCACGG"

#Spacer04/12
forward12="ACGTCATACGCCGACCAATCATCAG"
reverse12="CTGATGATTGGTCGGCGTATGACGT"

#Spacer05/13
forward13="TTTTCTGACCACTTGTGCGGGATTA"
reverse13="TAATCCCGCACAAGTGGTCAGAAAA"

#Spacer06/14
forward14="CGTCGTCATTTCCGGCTTCAATTTC"
reverse14="GAAATTGAAGCCGGAAATGACGACG"

#Spacer07/15
forward15="GAGGAGAGCGAGTACTCGGGGCTGC"
reverse15="GCAGCCCCGAGTACTCGCTCTCCTC"

#Spacer -/16
forward16="ACGTTAGGGCATGCAG"
reverse16="CTGCATGCCCTAACGT"

#Spacer -/17
forward17="TCTTGAGCAACGCCATCA"
reverse17="TGATGGCGTTGCTCAAGA"

#Spacer08/18
forward18="CGTGAAACCGCCCCCAGCCTCGCCG"
reverse18="CGGCGAGGCTGGGGGCGGTTTCACG"

#Spacer09/19
forward19="ACTCGGAATCCCATGTGCTGACAGC"
reverse19="GCTGTCAGCACATGGGATTCCGAGT"

#Spacer10/20
forward20="TCGACACCCGCTCTAGTTGACTTCC"
reverse20="GGAAGTCAACTAGAGCGGGTGTCGA"

#Spacer11/21
forward21="GTGAGCAACGGCGGCGGCAACCTGG"
reverse21="CCAGGTTGCCGCCGCCGTTGCTCAC"

#Spacer12/22
forward22="ATATCTGCTGCCCGCCCGGGGAGAT"
reverse22="ATCTCCCCGGGCGGGCAGCAGATAT"

#Spacer13/23
forward23="GACCATCATTGCCATTCCCTCTCCC"
reverse23="GGGAGAGGGAATGGCAATGATGGTC"

#Spacer14/24
forward24="GGTGTGATGCGGATGGTCGGCTCGG"
reverse24="CCGAGCCGACCATCCGCATCACACC"

#Spacer15/25
forward25="CTTGAATAACGCGCAGTGAATTTCG"
reverse25="CGAAATTCACTGCGCGTTATTCAAG"

#Spacer16/26
forward26="CGAGTTCCCGTCAGCGTCGTAAATC"
reverse26="GATTTACGACGCTGACGGGAACTCG"

#Spacer17/27
forward27="GCGCCGGCCCGCGCGGATGACTCCG"
reverse27="CGGAGTCATCCGCGCGGGCCGGCGC"

#Spacer18/28
forward28="CATGGACCCGGGCGAGCTGCAGATG"
reverse28="CATCTGCAGCTCGCCCGGGTCCATG"

#Spacer19/29
forward29="TAACTGGCTTGGCGCTGATCCTGGT"
reverse29="ACCAGGATCAGCGCCAAGCCAGTTA"

#Spacer20/30
forward30="TTGACCTCGCCAGGAGAGAAGATCA"
reverse30="TGATCTTCTCTCCTGGCGAGGTCAA"

#Spacer21/31
forward31="TCGATGTCGATGTCCCAATCGTCGA"
reverse31="TCGACGATTGGGACATCGACATCGA"

#Spacer22/32
forward32="ACCGCAGACGGCACGATTGAGACAA"
reverse32="TTGTCTCAATCGTGCCGTCTGCGGT"

#Spacer23/33
forward33="AGCATCGCTGATGCGGTCCAGCTCG"
reverse33="CGAGCTGGACCGCATCAGCGATGCT"

#Spacer24/34
forward34="CCGCCTGCTGGGTGAGACGTGCTCG"
reverse34="CGAGCACGTCTCACCCAGCAGGCGG"

#Spacer25/35
forward35="GATCAGCGACCACCGCACCCTGTCA"
reverse35="TGACAGGGTGCGGTGGTCGCTGATC"

#Spacer26/36
forward36="CTTCAGCACCACCATCATCCGGCGC"
reverse36="GCGCCGGATGATGGTGGTGCTGAAG"

#Spacer27/37
forward37="GGATTCGTGATCTCTTCCCGCGGAT"
reverse37="ATCCGCGGGAAGAGATCACGAATCC"

#Spacer28/38
forward38="TGCCCCGGCGTTTAGCGATCACAAC"
reverse38="GTTGTGATCGCTAAACGCCGGGGCA"

#Spacer29/39
forward39="AAATACAGGCTCCACGACACGACCA"
reverse39="TGGTCGTGTCGTGGAGCCTGTATTT"

#Spacer30/40
forward40="GGTTGCCCCGCGCCCTTTTCCAGCC"
reverse40="GGCTGGAAAAGGGCGCGGGGCAACC"

#Spacer31/41
forward41="TCAGACAGGTTCGCGTCGATCAAGT"
reverse41="ACTTGATCGACGCGAACCTGTCTGA"

#Spacer32/42
forward42="GACCAAATAGGTATCGGCGTGTTCA"
reverse42="TGAACACGCCGATACCTATTTGGTC"

#Spacer33/43
forward43="GACATGACGGCGGTGCCGCACTTGA"
reverse43="TCAAGTGCGGCACCGCCGTCATGTC"

#Spacer34/44
forward44="AAGTCACCTCGCCCACACCGTCGAA"
reverse44="TTCGACGGTGTGGGCGAGGTGACTT"

#Spacer -/45
forward45="AAGTTGGCGCTGGGG"
reverse45="CCCCAGCGCCAACTT"

#Spacer35/46
forward46="TCCGTACGCTCGAAACGCTTCCAAC"
reverse46="GTTGGAAGCGTTTCGAGCGTACGGA"

#Spacer36/47
forward47="CGAAATCCAGCACCACATCCGCAGC"
reverse47="GCTGCGGATGTGGTGCTGGATTTCG"

#Spacer -/48
forward48="AACCGTCCCACCTGG"
reverse48="CCAGGTGGGACGGTT"

#Spacer -/49
forward49="AACACTTTTTTTGAGCGTGG"
reverse49="CCACGCTCAAAAAAAGTGTT"

#Spacer -/50
forward50="CGGAAACGCAGCACC"
reverse50="GGTGCTGCGTTTCCG"

#Spacer37/51
forward51="CGCGAACTCGTCCACAGTCCCCCTT"
reverse51="AAGGGGGACTGTGGACGAGTTCGCG"

#Spacer38/52
forward52="CGTGGATGGCGGATGCGTTGTGCGC"
reverse52="GCGCACAACGCATCCGCCATCCACG"

#Spacer39/53
forward53="GACGATGGCCAGTAAATCGGCGTGG"
reverse53="CCACGCCGATTTACTGGCCATCGTC"

#Spacer -/54
forward54="CGATCATGAGAGTTGCG"
reverse54="CGCAACTCTCATGATCG"

#Spacer -/55
forward55="TTTTCGCTGTTGTGGTTCT"
reverse55="AGAACCACAACAGCGAAAA"

#Spacer -/56
forward56="AGCACCTCCCTTGACAA"
reverse56="TTGTCAAGGGAGGTGCT"

#Spacer -/57
forward57="TGCTGACTTCGCCTGTA"
reverse57="TACAGGCGAAGTCAGCA"

#Spacer -/58
forward58="CGAGCAGCGGCATA"
reverse58="TATGCCGCTGCTCG"

#Spacer -/59
forward59="GCATCCACTCGTCGC"
reverse59="GCGACGAGTGGATGC"

#Spacer -/60
forward60="TGGGTAATTGCGTCACGG"
reverse60="CCGTGACGCAATTACCCA"

#Spacer -/61
forward61="ACCATCCGACGCAGG"
reverse61="CCTGCGTCGGATGGT"

#Spacer40/62
forward62="CGCCATCTGTGCCTCATACAGGTCC"
reverse62="GGACCTGTATGAGGCACAGATGGCG"

#Spacer41/63
forward63="GGAGCTTTCCGGCTTCTATCAGGTA"
reverse63="TACCTGATAGAAGCCGGAAAGCTCC"

#Spacer42/64
forward64="ATGGTGGGACATGGACGAGCGCGAC"
reverse64="GTCGCGCTCGTCCATGTCCCACCAT"

#Spacer43/65
forward65="CGCAGAATCGCACCGGGTGCGGGAG"
reverse65="CTCCCGCACCCGGTGCGATTCTGCG"

#Spacer -/66
forward66="CCACGCTACTGCTCC"
reverse66="GGAGCAGTAGCGTGG"

#Spacer -/67
forward67="CACCGCCGATGACAG"
reverse67="CTGTCATCGGCGGTG"

#Spacer -/68
forward68="GTGTTTCGGCCGTGC"
reverse68="GCACGGCCGAAACAC"

#Spacer -/69
forward69="GTTGCATTCGTCGACTG"
reverse69="CAGTCGACGAATGCAAC"

#Spacer -/70
forward70="GGCGGCGCCGAGAA"
reverse70="TTCTCGGCGCCGCC"

#Spacer -/71
forward71="TTCCATGACTTGACGCC"
reverse71="GGCGTCAAGTCATGGAA"

#Spacer -/72
forward72="CGATGCGGCCACTAG"
reverse72="CTAGTGGCCGCATCG"

#Spacer -/73
forward73="GCTGACCCCATGGATG"
reverse73="CATCCATGGGGTCAGC"

#Spacer -/74
forward74="CAACAAGGTCTACGCGT"
reverse74="ACGCGTAGACCTTGTTG"

#Spacer -/75
forward75="GATCAGGCGAAGGCG"
reverse75="CGCCTTCGCCTGATC"

#Spacer -/76
forward76="ATTGCAGCGACGGGC"
reverse76="GCCCGTCGCTGCAAT"

#Spacer -/77
forward77="CAACGACGCTGTATT"
reverse77="AATACAGCGTCGTTG"

#Spacer -/78
forward78="AGCAGCATGGACGGTTT"
reverse78="AAACCGTCCATGCTGCT"

#Spacer -/79
forward79="GCGGATGTGGTGGTC"
reverse79="GACCACCACATCCGC"

#Spacer -/80
forward80="GTACATAGCGAGCTG"
reverse80="CAGCTCGCTATGTAC"

#Spacer -/81
forward81="GCCGCGGGTTTCGTT"
reverse81="AACGAAACCCGCGGC"

#Spacer -/82
forward82="GGGGCGTGTGTTCGTA"
reverse82="TACGAACACACGCCCC"

#Spacer -/83
forward83="CTGGTGTGCTTATGCCT"
reverse83="AGGCATAAGCACACCAG"

#Spacer -/84
forward84="CAAATGTTTGGACTGTGATC"
reverse84="GATCACAGTCCAAACATTTG"

#Spacer -/85
forward85="TTGTCGCGCGCCTTTTT"
reverse85="AAAAAGGCGCGCGACAA"

#Spacer -/86
forward86="GTTTCAGTTTTCTTGTCCC"
reverse86="GGGACAAGAAAACTGAAAC"

#Spacer -/87
forward87="CTGGTTGTTGCCCGG"
reverse87="CCGGGCAACAACCAG"

#Spacer -/88
forward88="TGTTCGGTGTTCTCTG"
reverse88="CAGAGAACACCGAACA"

#Spacer -/89
forward89="TCATGACGAGCCCGCA"
reverse89="TGCGGGCTCGTCATGA"

#Spacer -/90
forward90="ACACGGCCTGATCGGT"
reverse90="ACCGATCAGGCCGTGT"

#Spacer -/91
forward91="CGGATTGTCTGGCCC"
reverse91="GGGCCAGACAATCCG"

#Spacer -/92
forward92="TAAGCACGCGTCTGTCA"
reverse92="TGACAGACGCGTGCTTA"

#Spacer -/93
forward93="GACCACCGAATCACCAT"
reverse93="ATGGTGATTCGGTGGTC"

#Spacer -/94
forward94="TCTGGTAGTGGGCTTCT"
reverse94="AGAAGCCCACTACCAGA"


cat $forReads $revReads>seq

echo "**********Starting Spoligo Loop**********"

for i in {1..94}; do
	#echo "Currently on $i $(date +%F_%T)"
	forward=forward$i
	reverse=reverse$i
	sp=$i
	onemismatch
	#cursp=`cat $forReads $revReads | egrep -h -m 5 -c -f searchpattern`
	cursp=`egrep -h -m 5 -c -f searchpattern seq`
	extspacer="$extspacer	$cursp"
	if [ $( grep -wc $sp NUM) -gt 0 ]; then
	spacer="$spacer	$cursp"
	fi
	echo "Extended $extspacer, Original $spacer"
	rm searchpattern
	done



rm NUM
rm seq

echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43	spacer44	spacer 45	spacer46	spacer47	spacer48	spacer49	spacer50	spacer51	spacer52	spacer53	spacer54	spacer55	spacer56	spacer57	spacer58	spacer59	spacer60	spacer61	spacer62	spacer63	spacer64	spacer65	spacer66	spacer67 spacer68	spacer69	spacer70	spacer71	spacer72	spacer73	spacer74	spacer75	spacer76	spacer77	spacer78	spacer79	spacer80	spacer 81	spacer 82	spacer 83	spacer84	spacer85	spacer86	spacer87	spacer88	spacer89	spacer90	spacer91	spacer92	spacer93	spacer94" >> $n.spacer.txt
echo "$extspacer">$n.spacer.txt
echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43" >> $n.spacer.txt
echo "$spacer">>$n.spacer.txt

mybinaries=`echo $spacer | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
#`echo $spacer | awk 'NR==2 {for(i=1;i<=NF;i++) if ($i >= 5) print 1; else print 0}' | tr -cd "[:print:]" | fold -w3`
myextbinaries=`echo $extspacer | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
echo "Binaries: $mybinaries, $myextbinaries"

exspoligooctal=`echo "ibase=2;obase=8; $mybinaries" | bc`>$n.octalcode.txt
spoligooctal=`echo "ibase=2;obase=8; $myextbinaries" | bc`>>$n.octalcode.txt
echo "Octals: $exspoligooctal, $spoligooctal"

#for i in $mybinaries; do 
#if [ $i == 000 ]
#then
#	echo "0" >> $n.octalcode.txt
#elif [ $i == 001 ]
#then
#	echo "1" >> $n.octalcode.txt
#elif [ $i == 010 ]
#	then
#	echo "2" >> $n.octalcode.txt
#elif [ $i == 011 ]
#	then
#	echo "3" >> $n.octalcode.txt
#elif [ $i == 100 ]
#	then
#	echo "4" >> $n.octalcode.txt
#elif [ $i == 101 ]
#	then
#	echo "5" >> $n.octalcode.txt
#elif [ $i == 110 ]
#	then
#	echo "6" >> $n.octalcode.txt	
#elif [ $i == 111 ]
#	then	
#	echo "7" >> $n.octalcode.txt
#elif [ $i == 0 ]
#	then	
#	echo "0" >> $n.octalcode.txt
#elif [ $i == 1 ]
#	then
#	echo "1" >> $n.octalcode.txt
#else
#	echo "***Error***" >> $n.octalcode.txt
#fi
#done

tr -d '\n' < $n.octalcode.txt | awk -v x="$n" 'BEGIN{OFS="\t"}{print x, $0}' > spoligo.txt
WGSpoligo=`cat $n.octalcode.txt | tr -cd "[:print:]"`

# Add infor to spoligoCheck.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck.txt
echo "WGExtSpoligo:	$exspoligooctal" >> /scratch/report/spoligoCheck.txt

#Make FileMaker file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "$n" "$WGSpoligo" >> "/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

# Add infor to spoligoCheck_all.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck_all.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck_all.txt
echo "WGExtSpoligo:	$exspoligooctal" >> /scratch/report/spoligoCheck_all.txt

# move back a directory to main sample folder
cd ..

#
#  Created by Stuber, Tod P - APHIS on 03/07/2013.
#  Update by Hicks, Jessica A - APHIS on 01/23/2017.
