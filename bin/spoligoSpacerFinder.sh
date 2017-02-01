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
#Original sequences are 25bp probe, but new sequences had shorter probes.  Used additional spacer sequence to extend to 25 bp

#Spacer -/01
forward1="GCACTGCAACCCGGAATTCTTGCAC"
reverse1="GTGCAAGAATTCCGGGTTGCAGTGC"

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
forward5="ATGGCACGGCAGGCGTGGCTAGGGG"
reverse5="CCCCTAGCCACGCCTGCCGTGCCAT"

#Spacer -/6
forward6="ATGTGCGCCGTCGCCGTAAGTGCCC"
reverse6="GGGCACTTACGGCGACGGCGCACAT"

#Spacer -/7
forward7="AATTCGTTGACCACGAATTTTCAGA"
reverse7="TCTGAAAATTCGTGGTCAACGAATT"

#Spacer -/8
forward8="ACCGCTGGCGCGCATCATTCATCGA"
reverse8="TCGATGAATGATGCGCGCCAGCGGT"

#Spacer -/9
forward9="CCATATCGGGGACGGCGACGCTGCG"
reverse9="CGCAGCGTCGCCGTCCCCGATATGG"

#Spacer -/10
forward10="TACACCACGCGTCGTGCCATCAGTC"
reverse10="GACTGATGGCACGACGCGTGGTGTA"

#Spacer -/11
forward11="CCGTGCACATGCCGTGGCTCAGGGG"
reverse11="CCCCTGAGCCACGGCATGTGCACGG"

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
forward16="ACGACGTTAGGGCATGCAGCATGCC"
reverse16="GGCATGCTGCATGCCCTAACGTCGT"

#Spacer -/17
forward17="TGCTCTTGAGCAACGCCATCATCCG"
reverse17="CGGATGATGGCGTTGCTCAAGAGCA"

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
forward45="AAGTTGGCGCTGGGGTCTGAGTCAA"
reverse45="TTGACTCAGACCCCAGCGCCAACTT"

#Spacer35/46
forward46="TCCGTACGCTCGAAACGCTTCCAAC"
reverse46="GTTGGAAGCGTTTCGAGCGTACGGA"

#Spacer36/47
forward47="CGAAATCCAGCACCACATCCGCAGC"
reverse47="GCTGCGGATGTGGTGCTGGATTTCG"

#Spacer -/48
forward48="GCGAGGAACCGTCCCACCTGGGCCT"
reverse48="AGGCCCAGGTGGGACGGTTCCTCGC"

#Spacer -/49
forward49="TCAATAACACTTTTTTTGAGCGTGG"
reverse49="CCACGCTCAAAAAAAGTGTTATTGA"

#Spacer -/50
forward50="ACGGAAACGCAGCACCAGCCTGACA"
reverse50="TGTCAGGCTGGTGCTGCGTTTCCGT"

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
forward54="CATCGATCATGAGAGTTGCGTTGAT"
reverse54="ATCAACGCAACTCTCATGATCGATG"

#Spacer -/55
forward55="GATTTTCGCTGTTGTGGTTCTCATT"
reverse55="AATGAGAACCACAACAGCGAAAATC"

#Spacer -/56
forward56="GCACACCAGCACCTCCCTTGACAAT"
reverse56="ATTGTCAAGGGAGGTGCTGGTGTGC"

#Spacer -/57
forward57="CCTAAGGGTGCTGACTTCGCCTGTA"
reverse57="TACAGGCGAAGTCAGCACCCTTAGG"

#Spacer -/58
forward58="CCGACGACCGAGCAGCGGCATAGA"
reverse58="TCTATGCCGCTGCTCGGTCGTCGG"

#Spacer -/59
forward59="TTGCATCCACTCGTCGCCGACACGG"
reverse59="CCGTGTCGGCGACGAGTGGATGCAA"

#Spacer -/60
forward60="TGGGTAATTGCGTCACGGCGCGCCTG"
reverse60="CAGGCGCGCCGTGACGCAATTACCCA"

#Spacer -/61
forward61="ACCATCCGACGCAGGCACCGAAGTC"
reverse61="GACTTCGGTGCCTGCGTCGGATGGT"

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
forward66="CACCACAGCCACGCTACTGCTCCAT"
reverse66="ATGGAGCAGTAGCGTGGCTGTGGTG"

#Spacer -/67
forward67="ACACCGCCGATGACAGCTATGTCCG"
reverse67="CGGACATAGCTGTCATCGGCGGTGT"

#Spacer -/68
forward68="CGCGCGGTGTTTCGGCCGTGCCCGA"
reverse68="TCGGGCACGGCCGAAACACCGCGCG"

#Spacer -/69
forward69="GTTGCATTCGTCGACTGCGTGGTAT"
reverse69="ATACCACGCAGTCGACGAATGCAAC"

#Spacer -/70
forward70="GTAGCGGCCCCGGCGGCGCCGAGAA"
reverse70="TTCTCGGCGCCGCCGGGGCCGCTAC"

#Spacer -/71
forward71="TGGTGATCTTCCATGACTTGACGCC"
reverse71="GGCGTCAAGTCATGGAAGATCACCA"

#Spacer -/72
forward72="GCGGTGCTCGATGCGGCCACTAGGC"
reverse72="GCCTAGTGGCCGCATCGAGCACCGC"

#Spacer -/73
forward73="TCGGTGCTGACCCCATGGATGCGAA"
reverse73="TTCGCATCCATGGGGTCAGCACCGA"

#Spacer -/74
forward74="CAACAAGGTCTACGCGTCGAGGTCC"
reverse74="GGACCTCGACGCGTAGACCTTGTTG"

#Spacer -/75
forward75="ATTACGCCTGATCAGGCGAAGGCGA"
reverse75="TCGCCTTCGCCTGATCAGGCGTAAT"

#Spacer -/76
forward76="TTCAGTAAATTGCAGCGACGGGCGA"
reverse76="TCGCCCGTCGCTGCAATTTACTGAA"

#Spacer -/77
forward77="CTTCAACGACGCTGTATTGGGCCAT"
reverse77="ATGGCCCAATACAGCGTCGTTGAAG"

#Spacer -/78
forward78="AGCAGCATGGACGGTTTCGCCTGTA"
reverse78="TACAGGCGAAACCGTCCATGCTGCT"

#Spacer -/79
forward79="GTTGCGGATGTGGTGGTCGCGTAGC"
reverse79="GCTACGCGACCACCACATCCGCAAC"

#Spacer -/80
forward80="TTGGCGTACATAGCGAGCTGTGCGG"
reverse80="CCGCACAGCTCGCTATGTACGCCAA"

#Spacer -/81
forward81="TTGTGCCGCCGCGGGTTTCGTTCAC"
reverse81="GTGAACGAAACCCGCGGCGGCACAA"

#Spacer -/82
forward82="GGGGCGTGTGTTCGTAGTCGCCTAA"
reverse82="TTAGGCGACTACGAACACACGCCCC"

#Spacer -/83
forward83="GTGCTGGTGTGCTTATGCCTAACAG"
reverse83="CTGTTAGGCATAAGCACACCAGCAC"

#Spacer -/84
forward84="CAAATGTTTGGACTGTGATCAATTC"
reverse84="GAATTGATCACAGTCCAAACATTTG"

#Spacer -/85
forward85="TTGTCGCGCGCCTTTTTCCAGCCGA"
reverse85="TCGGCTGGAAAAAGGCGCGCGACAA"

#Spacer -/86
forward86="GCGTTTCAGTTTTCTTGTCCCAGTG"
reverse86="CACTGGGACAAGAAAACTGAAACGC"

#Spacer -/87
forward87="ACTGGTTGTTGCCCGGCGACGGCGG"
reverse87="CCGCCGTCGCCGGGCAACAACCAGT"

#Spacer -/88
forward88="AAGTGGTGTTCGGTGTTCTCTGTAC"
reverse88="GTACAGAGAACACCGAACACCACTT"

#Spacer -/89
forward89="CGATCCGGTCATGACGAGCCCGCAG"
reverse89="CTGCGGGCTCGTCATGACCGGATCG"

#Spacer -/90
forward90="ATCACGACACGGCCTGATCGGTGTC"
reverse90="GACACCGATCAGGCCGTGTCGTGAT"

#Spacer -/91
forward91="GCGTCGGCCGGATTGTCTGGCCCAC"
reverse91="GTGGGCCAGACAATCCGGCCGACGC"

#Spacer -/92
forward92="CGTCGGCTAAGCACGCGTCTGTCAA"
reverse92="TTGACAGACGCGTGCTTAGCCGACG"

#Spacer -/93
forward93="GGTGAGGACCACCGAATCACCATCA"
reverse93="TGATGGTGATTCGGTGGTCCTCACC"

#Spacer -/94
forward94="TCTGGTAGTGGGCTTCTGCCGGTGC"
reverse94="GCACCGGCAGAAGCCCACTACCAGA"


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
	#echo "Extended $extspacer, Original $spacer"
	rm searchpattern
	done



rm NUM
rm seq

echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43	spacer44	spacer 45	spacer46	spacer47	spacer48	spacer49	spacer50	spacer51	spacer52	spacer53	spacer54	spacer55	spacer56	spacer57	spacer58	spacer59	spacer60	spacer61	spacer62	spacer63	spacer64	spacer65	spacer66	spacer67 spacer68	spacer69	spacer70	spacer71	spacer72	spacer73	spacer74	spacer75	spacer76	spacer77	spacer78	spacer79	spacer80	spacer 81	spacer 82	spacer 83	spacer84	spacer85	spacer86	spacer87	spacer88	spacer89	spacer90	spacer91	spacer92	spacer93	spacer94" > $n.spacer.txt
echo "$extspacer" >> $n.spacer.txt
echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43" >> $n.spacer.txt
echo "$spacer" >> $n.spacer.txt


mybinaries=`echo $spacer | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
myextbinaries=`echo $extspacer | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]"`
#echo "Binaries: $mybinaries, $myextbinaries"

#remove last digit to make divisible by 3 for the bc function to convert to octal, spoligo reads left to right, but bc converts right to left
spbinary=`echo $mybinaries | awk '{print substr($1,1,42)}'`
spbinarylast=`echo $mybinaries | awk '{print substr($1,43,1)}'`
exspbinary=`echo $myextbinaries | awk '{print substr($1,1,93)}'`
exspbinarylast=`echo $myextbinaries | awk '{print substr($1,94,1)}'`

#convert to octal and add last digit back
spoligooctal=`echo "ibase=2;obase=8; $spbinary" | bc`
WGSpoligo=`echo "$spoligooctal$spbinarylast"`
echo "$WGSpoligo">$n.octalcode.txt
exspoligooctal=`echo "ibase=2;obase=8; $exspbinary" | bc`
WGExSpoligo=`echo "$exspoligooctal$exspbinarylast"`
echo "$WGExSpoligo">>$n.octalcode.txt
#echo "Octals: $exspoligooctal, $spoligooctal"


echo "$n  $WGSpoligo">spoligo.txt
echo "$WGExSpoligo">>spoligo.txt
# Add infor to spoligoCheck.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck.txt
echo "WGExSpoligo:	$WGExSpoligo" >> /scratch/report/spoligoCheck.txt

#Make FileMaker file
dateFile=`date "+%Y%m%d"`
printf "%s\t%s\n" "$n" "$WGSpoligo" >> "/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/newFiles/${dateFile}_FileMakerSpoligoImport.txt"

# Add infor to spoligoCheck_all.txt
echo "<----- $n ----->" >> /scratch/report/spoligoCheck_all.txt
echo "WGSpoligo:	$WGSpoligo" >> /scratch/report/spoligoCheck_all.txt
echo "WGExtSpoligo:	$WGExSpoligo" >> /scratch/report/spoligoCheck_all.txt

# move back a directory to main sample folder
cd ..

#
#  Created by Stuber, Tod P - APHIS on 03/07/2013.
#  Update by Hicks, Jessica A - APHIS on 01/23/2017.
