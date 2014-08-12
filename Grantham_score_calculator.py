## This program takes input a SnpEff annotated vcf file and calcuates Grantham substitution scores for the missense variants.


import math,sys,re,fileinput,os

'''
Granthan score matrix 
		Arg	Leu	Pro	Thr	Ala	Val	Gly	Ile	Phe	Tyr	Cys	His	Gln	Asn	Lys	Asp	Glu	Met	Trp		
		R	L	P	T	A	V	G	I	F	Y	C	H	Q	N	K	D	E	M	W		
Ser	S	110	145	74	58	99	124	56	142	155	144	112	89	68	46	121	65	80	135	177	S	Ser
Arg	R	0	102	103	71	112	96	125	97	97	77	180	29	43	86	26	96	54	91	101	R	Arg
Leu	L	0	0	98	92	96	32	138	5	22	36	198	99	113	153	107	172	138	15	61	L	Leu
Pro	P	0	0	0	38	27	68	42	95	114	110	169	77	76	91	103	108	93	87	147	P	Pro
Thr	T	0	0	0	0	58	69	59	89	103	92	149	47	42	65	78	85	65	81	128	T	Thr
Ala	A	0	0	0	0	0	64	60	94	113	112	195	86	91	111	106	126	107	84	148	A	Ala
Val	V	0	0	0	0	0	0	109	29	50	55	192	84	96	133	97	152	121	21	88	V	Val
Gly	G	0	0	0	0	0	0	0	135	153	147	159	98	87	80	127	94	98	127	184	G	Gly
Ile	I	0	0	0	0	0	0	0	0	21	33	198	94	109	149	102	168	134	10	61	I	Ile
Phe	F	0	0	0	0	0	0	0	0	0	22	205	100	116	158	102	177	140	28	40	F	Phe
Tyr	Y	0	0	0	0	0	0	0	0	0	0	194	83	99	143	85	160	122	36	37	Y	Tyr
Cys	C	0	0	0	0	0	0	0	0	0	0	0	174	154	139	202	154	170	196	215	C	Cys
His	H	0	0	0	0	0	0	0	0	0	0	0	0	24	68	32	81	40	87	115	H	His
Gln	Q	0	0	0	0	0	0	0	0	0	0	0	0	0	46	53	61	29	101	130	Q	Gln
Asn	N	0	0	0	0	0	0	0	0	0	0	0	0	0	0	94	23	42	142	174	N	Asn
Lys	K	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	101	56	95	110	K	Lys
Asp	D	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	160	181	D	Asp
Glu	E	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	126	152	E	Glu
Met	M	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	67	M	Met

'''

def average(x):
	assert len(x) > 0
	return float(sum(x))/float(len(x))

Argument = []
Argument = sys.argv[1:] 



if (len(Argument)) < 2:	
	print "Usage: Input_variant_file Output_file" 
	sys.exit()

Filepath = Argument[0]
output = open(Argument[1],"w")

def Grantham_score_transcript(variant_list):
	 
	Sum = 0
	variants = []
	variants = variant_list

	Grantham_matrix = {

'S':	{'R':110,	'L':145,	'P':74,	'T':58,	'A':99,	'V':124,	'G':56,	'I':142,	'F':155,	'Y':144,	'C':112,	'H':89,	'Q':68,	'N':46,	'K':121,	'D':65,	'E':80,	'M':135,	'W':177},
'R':	{'R':0,	'L':102,	'P':103,	'T':71,	'A':112,	'V':96,	'G':125,	'I':97,	'F':97,	'Y':77,	'C':180,	'H':29,	'Q':43,	'N':86,	'K':26,	'D':96,	'E':54,	'M':91,	'W':101, 'S':0},
'L':	{'R':0,	'L':0,	'P':98,	'T':92,	'A':96,	'V':32,	'G':138,	'I':5,	'F':22,	'Y':36,	'C':198,	'H':99,	'Q':113,	'N':153,	'K':107,	'D':172,	'E':138,	'M':15,	'W':61, 'S':0},
'P':	{'R':0,	'L':0,	'P':0,	'T':38,	'A':27,	'V':68,	'G':42,	'I':95,	'F':114,	'Y':110,	'C':169,	'H':77,	'Q':76,	'N':91,	'K':103,	'D':108,	'E':93,	'M':87,	'W':147, 'S':0},
'T':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':58,	'V':69,	'G':59,	'I':89,	'F':103,	'Y':92,	'C':149,	'H':47,	'Q':42,	'N':65,	'K':78,	'D':85,	'E':65,	'M':81,	'W':128, 'S':0},
'A':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':64,	'G':60,	'I':94,	'F':113,	'Y':112,	'C':195,	'H':86,	'Q':91,	'N':111,	'K':106,	'D':126,	'E':107,	'M':84,	'W':148, 'S':0},
'V':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':109,	'I':29,	'F':50,	'Y':55,	'C':192,	'H':84,	'Q':96,	'N':133,	'K':97,	'D':152,	'E':121,	'M':21,	'W':88, 'S':0},
'G':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':135,	'F':153,	'Y':147,	'C':159,	'H':98,	'Q':87,	'N':80,	'K':127,	'D':94,	'E':98,	'M':127,	'W':184, 'S':0},
'I':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':21,	'Y':33,	'C':198,	'H':94,	'Q':109,	'N':149,	'K':102,	'D':168,	'E':134,	'M':10,	'W':61, 'S':0},
'F':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':22,	'C':205,	'H':100,	'Q':116,	'N':158,	'K':102,	'D':177,	'E':140,	'M':28,	'W':40, 'S':0},
'Y':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':194,	'H':83,	'Q':99,	'N':143,	'K':85,	'D':160,	'E':122,	'M':36,	'W':37, 'S':0},
'C':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':174,	'Q':154,	'N':139,	'K':202,	'D':154,	'E':170,	'M':196,	'W':215, 'S':0},
'H':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':24,	'N':68,	'K':32,	'D':81,	'E':40,	'M':87,	'W':115, 'S':0},
'Q':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':46,	'K':53,	'D':61,	'E':29,	'M':101,	'W':130, 'S':0},
'N':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':0,	'K':94,	'D':23,	'E':42,	'M':142,	'W':174, 'S':0},
'K':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':0,	'K':0,	'D':101,	'E':56,	'M':95,	'W':110, 'S':0},
'D':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':0,	'K':0,	'D':0,	'E':45,	'M':160,	'W':181, 'S':0},
'E':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':0,	'K':0,	'D':0,	'E':0,	'M':126,	'W':152, 'S':0},
'M':	{'R':0,	'L':0,	'P':0,	'T':0,	'A':0,	'V':0,	'G':0,	'I':0,	'F':0,	'Y':0,	'C':0,	'H':0,	'Q':0,	'N':0,	'K':0,	'D':0,	'E':0,	'M':0,	'W':67, 'S':0}, 
'W':    {'R':0, 'L':0,  'P':0,  'T':0,  'A':0,  'V':0,  'G':0,  'I':0,  'F':0,  'Y':0,  'C':0,  'H':0,  'Q':0,  'N':0,  'K':0,  'D':0,  'E':0,  'M':0,  'W':0 , 'S':0}}

	
	for variant in variants:
		Ref_aa = ""
		Alt_aa = ""

		Ref_aa = str(variant[0])
		Alt_aa = str(variant[-1])

		#print Ref_aa
		#print Alt_aa

		if Ref_aa == Alt_aa:
			return Sum
		else:
			if int(Grantham_matrix[Ref_aa][Alt_aa]) != 0:
				Sum = Sum+int(Grantham_matrix[Ref_aa][Alt_aa])
			else:
				Sum = Sum+int(Grantham_matrix[Alt_aa][Ref_aa])

	if Sum > 150:	
		return str(Sum)+"\t"+"Radical"
	elif Sum <= 150 and Sum > 100:
		return str(Sum)+"\t"+"Moderately Radical"
	elif Sum <= 100 and Sum > 50:
		return str(Sum)+"\t"+"Moderately Conservative"
	else:
		return str(Sum)+"\t"+"Conservative"
	
Transcripts_total = {}
Info = {}

for line in fileinput.input([Filepath]):
	rowlist = []
	rowlist = (line.rstrip("\n")).split('\t')
	
	if line.startswith("#"):
		continue
	else:
		temp1 = []
		
		temp1 = ((rowlist[7].split(";")[-1]).lstrip("EFF=")).split(",")
		for effect in temp1:
			if effect.startswith("NON_SYNONYMOUS_CODING"):
				temp2 = []
				temp2 = (effect.lstrip("NON_SYNONYMOUS_CODING\(MODERATE|MISSENSE|").rstrip("\)")).split("|")
				#print temp2
	  
			if temp2[-2] not in Transcripts_total:
				Transcripts_total[temp2[-2]] = []
				Transcripts_total[temp2[-2]].append(temp2[1])
			else:
				Transcripts_total[temp2[-2]].append(temp2[1])

for line1 in fileinput.input([Filepath]):
        rowlist1 = []
        rowlist1 = (line1.rstrip("\n")).split('\t')

        if line1.startswith("#"):
                #output.write(str(line1))
                continue
        else:
                temp3 = []
                temp4 = []

                temp3 = ((rowlist1[7].split(";")[-1]).lstrip("EFF=")).split(",")
                for effect1 in temp3:
                        if effect1.startswith("NON_SYNONYMOUS_CODING"):
                                temp4 = (effect1.lstrip("NON_SYNONYMOUS_CODING\(MODERATE|MISSENSE|").rstrip("\)")).split("|")
				function_pass = []
				function_pass.append(temp4[1])

				#print line1

				output.write(str("\t".join(rowlist1[:5]))+"\t"+str(rowlist1[7].split(";")[-2])+"\t"+str("\t".join(temp4))+"\t"+str(Grantham_score_transcript(function_pass))+"\t"+str(Grantham_score_transcript(Transcripts_total[temp4[-2]])))

		output.write("\n")

output.close()


