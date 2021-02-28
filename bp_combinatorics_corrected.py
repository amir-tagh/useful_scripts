#!/usr/bin/env python3

import numpy as np

beta    = 1./0.59616108
n_samps = 8000000

##Data from paper by barry honig
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1236383/?page=5

##second set of data from Maxim Kamenetskii
##JMB, 2004, 342, 775-785

##Evolutionary stages of amino acids (2ndary source)
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3293468/
bases = ['A', 'T', 'G', 'C']

##chou & fasman 1978
##also arxiv... recent...
prop_alpha = {
'A' : 1.46, 'C' : 0.75,
'D' : 0.83, 'E' : 1.37,
'F' : 0.97, 'G' : 0.43,
'H' : 0.90,
'I' : 1.07, 'K' : 1.14,
'L' : 1.36, 'M' : 1.25,
'N' : 0.72, 'P' : 0.40,
'Q' : 1.34, 'R' : 1.24,
'S' : 0.76, 'T' : 0.76,
'V' : 0.94, 'W' : 1.89,
'Y' : 0.95, 
'X' : 0.
}
prop_beta = {
'A' : 0.78, 'C' : 1.26,
'D' : 0.51, 'E' : 0.69,
'F' : 1.49, 'G' : 0.69,
'H' : 1.05,
'I' : 1.65, 'K' : 0.78,
'L' : 1.13, 'M' : 1.10,
'N' : 0.64, 'P' : 0.32,
'Q' : 0.82, 'R' : 0.88,
'S' : 0.88, 'T' : 1.21,
'V' : 1.89, 'W' : 1.39,
'Y' : 1.47, 
'X' : 0.
}

##from genscript.com
codon_freq_e_coli = {
'TTT' : 0.0221, 'TAT' : 0.0175, 
'TTC' : 0.0160, 'TAC' : 0.0122,
'TTA' : 0.0143, 'TAA' : 0.0020,
'TTG' : 0.0130, 'TAG' : 0.0003,
'CTT' : 0.0119, 'CAT' : 0.0125,
'CTC' : 0.0102, 'CAC' : 0.0093,
'CTA' : 0.0042, 'CAA' : 0.0146,
'CTG' : 0.0484, 'CAG' : 0.0284,
'ATT' : 0.0298, 'AAT' : 0.0206,
'ATC' : 0.0237, 'AAC' : 0.0214,
'ATA' : 0.0068, 'AAA' : 0.0353,
'ATG' : 0.0264, 'AAG' : 0.0124,
'GTT' : 0.0198, 'GAT' : 0.0327,
'GTC' : 0.0143, 'GAC' : 0.0192,
'GTA' : 0.0116, 'GAA' : 0.0391,
'GTG' : 0.0244, 'GAG' : 0.0187,
'TCT' : 0.0104, 'TGT' : 0.0052,
'TCC' : 0.0091, 'TGC' : 0.0061,
'TCA' : 0.0089, 'TGA' : 0.0010,
'TCG' : 0.0085, 'TGG' : 0.0139,
'CCT' : 0.0075, 'CGT' : 0.0200,
'CCC' : 0.0054, 'CGC' : 0.0197,
'CCA' : 0.0086, 'CGA' : 0.0038,
'CCG' : 0.0209, 'CGG' : 0.0059,
'ACT' : 0.0103, 'AGT' : 0.0099,
'ACC' : 0.0220, 'AGC' : 0.0152,
'ACA' : 0.0093, 'AGA' : 0.0036,
'ACG' : 0.0137, 'AGG' : 0.0021,
'GCT' : 0.0171, 'GGT' : 0.0255,
'GCC' : 0.0242, 'GGC' : 0.0271,
'GCA' : 0.0212, 'GGA' : 0.0095,
'GCG' : 0.0301, 'GGG' : 0.0113,
}

t = 0.
keyList_ec = []
wList_ec   = []
for c in codon_freq_e_coli.keys():
    t += codon_freq_e_coli[c]
    keyList_ec.append(c)
    wList_ec.append(t)
print("total codon freq: %e" % t)




wc    = {
'C' : 'G',
'G' : 'C',
'A' : 'T',
'T' : 'A',
}

honig_data = {
'AA' : -6.53,
'AT' : -5.32,
'AG' : -7.03,
'AC' : -5.00,
'TA' : -5.25,
'TT' : -4.92,
'TG' : -5.25,
'TC' : -4.78,
'GA' : -7.41,
'GT' : -5.93,
'GG' : -7.79,
'GC' : -5.76,
'CA' : -4.92,
'CT' : -4.64,
'CG' : -5.07,
'CC' : -4.36,
}

kamanet_data ={
'AA' : -1.11,
'AT' : -1.34,
'AG' : -1.06,
'AC' : -1.81,
'TA' : -0.19,
'TT' : -1.11,
'TG' : -0.55,
'TC' : -1.43,
'GA' : -1.43,
'GT' : -1.81,
'GG' : -1.44,
'GC' : -2.17,
'CA' : -0.55,
'CT' : -1.06,
'CG' : -0.91,
'CC' : -1.44,
}
 
for k in honig_data:
    u          = honig_data[k]
    complement = wc[k[1]]+wc[k[0]]
    u         += honig_data[complement]
    print (k+"-"+complement+" "+str(u))


code = {
'GTA' : 'V',
'GTT' : 'V',
'GTG' : 'V',
'GTC' : 'V',
'GCA' : 'A',
'GCT' : 'A',
'GCG' : 'A',
'GCC' : 'A',
'GAA' : 'E',
'GAT' : 'D',
'GAG' : 'E',
'GAC' : 'D',
'GGA' : 'G',
'GGT' : 'G',
'GGG' : 'G',
'GGC' : 'G',
'ATA' : 'I',
'ATT' : 'I',
'ATG' : 'M',
'ATC' : 'I',
'ACA' : 'T',
'ACT' : 'T',
'ACG' : 'T',
'ACC' : 'T',
'AAA' : 'K',
'AAT' : 'N',
'AAG' : 'K',
'AAC' : 'N',
'AGA' : 'R',
'AGT' : 'S',
'AGG' : 'R',
'AGC' : 'S',
'CTA' : 'L',
'CTT' : 'L',
'CTG' : 'L',
'CTC' : 'L',
'CCA' : 'P',
'CCT' : 'P',
'CCG' : 'P',
'CCC' : 'P',
'CAA' : 'Q',
'CAT' : 'H',
'CAG' : 'Q',
'CAC' : 'H',
'CGA' : 'R',
'CGT' : 'R',
'CGG' : 'R',
'CGC' : 'R',
'TTA' : 'L',
'TTT' : 'F',
'TTG' : 'L',
'TTC' : 'F',
'TCA' : 'S',
'TCT' : 'S',
'TCG' : 'S',
'TCC' : 'S',
'TAA' : 'X',
'TAT' : 'Y',
'TAG' : 'X',
'TAC' : 'Y',
'TGA' : 'X',
'TGT' : 'C',
'TGG' : 'Y',
'TGC' : 'C',
}
reverse_code = {}
for codon in code.keys():
    aa = code[codon]
    if aa not in reverse_code:
        reverse_code[aa] = []
    reverse_code[aa].append(codon)


phase_one_aa = {
'L' : 1,
'S' : 1,
'P' : 1,
'I' : 1,
'T' : 1,
'V' : 1,
'A' : 1,
'D' : 1,
'E' : 1,
'G' : 1,
}
phase_zero_aa = {
'V' : 1,
'A' : 1,
'D' : 1,
'G' : 1,
}
p1_aa_list = []
for k in phase_one_aa.keys():
    p1_aa_list.append(k)


def get_random_codon( codset='ecoli' ):
    if codset == 'ecoli':
        w = np.random.random((1))
        i = 0
        while w > wList_ec[i]:
            i+=1;
            if i >= len(wList_ec):
                print("error picking a codon")
                exit(1)
    if codset == 'any_phase1':
        while True:
            i     = np.random.randint(64)
            codon = keyList_ec[i]
            aa    = code[codon]
            if aa in phase_one_aa:
                return codon
    if codset == 'phase1_GxC':
        i     = np.random.randint(len(p1_aa_list))
        aa    = p1_aa_list[i]
        codon = reverse_code[aa][0]
        cVar  = codon[:2]+'C'
        if cVar in reverse_code[aa]:
                codon = cVar
        gVar = 'G'+codon[1:]
        if gVar in reverse_code[aa]:
                codon = gVar
        return codon
    else:
        i = np.random.randint(64)
        return keyList_ec[i]

def get_utrip_inRandomSeq( b1, b2, b3 ):
    
    trip_prev = get_random_codon(codset='phase1_GxC')
    trip_next = get_random_codon(codset='phase1_GxC')

    ##looking at two strands of 3 codons here.
    ##: 8 possible breakage sites.
    s1=trip_prev+b1+b2+b3+trip_next
    s2=''
    for i in range(len(s1)):
        s2 += wc[s1[len(s1)-i-1]]

    ##collect energies to partition at any two sites:
    UtList = []
    UtBest = -9e9
    iBest  = 0

    partition = 0.
    for i in range(2,5):
        Ut1  = honig_data[s1[i:i+2]] + honig_data[s2[9-i-2:9-i]]
        for j in range(i+1,6):

            Ut2  = honig_data[s1[j:j+2]] + honig_data[s2[9-j-2:9-j]]
            U = Ut1+Ut2
            UtList.append( (i,j,U) )
            partition += np.exp(U * beta)
            if U > UtBest and not (i==2 and j==5):
                UtBest = U
                iBest  = len(UtList)-1
            elif (i==2 and j==5):
                cost = -1. * U

    ##return tiplet partition cost and 
    ##cost relative to most favourable alternative
    relCost = cost + UtBest;

    #print("partition function: ")
    #print(UtList)
    #print("%e : %e kT" % (partition, -1*np.log(partition) ))

    deltaG  = -1*np.log( np.exp(-1*cost*beta)/partition )

    #print("relative cost:"+str(relCost))

    return cost, relCost, deltaG

def get_utrip( b1, b2, b3 ):
    c3 = wc[b1]
    c2 = wc[b2]
    c1 = wc[b3]

    ##energy to break  a repeat sequence up into 
    ##triplets
    Ut  = honig_data[b1+b2] + honig_data[b2+b3]
    #Ut -= honig_data[b3+b1]             

    ##complement sequence must break at the same point
    Ut += honig_data[c1+c2] + honig_data[c2+c3]
    #Ut -= honig_data[c3+c1]   

    return Ut, c1, c2, c3



Utrip  = []
Utrip_random = []
G_random     = []
trips  = []
Ucomp  = []

AAlist = []
done = {}

for b1 in bases:
    for b2 in bases:
        for b3 in bases:

            if b1+b2+b3 in done:
                continue

            Ut, c1, c2, c3 = get_utrip( b1, b2, b3 )

            if c1+c2+c3 in done:
                continue

            done[b1+b2+b3] = 1
            done[c1+c2+c3] = 1

            Utrip.append( Ut )          

            ##try sampling N random neighbouring codons
            UtR = 0.
            UG  = 0.
            for i in range(n_samps):
                U, rel, G = get_utrip_inRandomSeq( b1, b2, b3 )
                UG   += G
                UtR  += rel
            UtR /= n_samps
            UG  /= n_samps
            Utrip_random.append( UtR )
            G_random.append(UG)

            ##competitor energies
            Uc1 = get_utrip( b2, b3, b1 )[0]
            Uc2 = get_utrip( b3, b1, b2 )[0]

            Uc  = min(Uc1, Uc2)
            Ucomp.append(Uc)

            if code[b1+b2+b3] in phase_one_aa:
                if code[b1+b2+b3] in phase_zero_aa:
                    phase1 = "** "
                else:
                    phase1 = "*  "
            else:
                phase1 = ".  "
            if code[c1+c2+c3] in phase_one_aa:
                if code[c1+c2+c3] in phase_zero_aa:
                    phase2 = "** "
                else:
                    phase2 = "*  "
            else:
                phase2 = ".  "
                

            trips.append(b1+b2+b3+" . "+c3+c2+c1+\
                         "  "+code[b1+b2+b3]+phase1+\
                         "  "+code[c1+c2+c3]+phase2)
            AAlist.append( [code[b1+b2+b3], code[c1+c2+c3]] )


Utrip        = np.array(Utrip)
Utrip_random = np.array(Utrip_random)
G_random     = np.array(G_random)
Ucomp = np.array(Ucomp)
Udelt = Utrip - Ucomp
index = np.argsort(G_random)


##find the best available triplet for each AA
best_U  = {}
best_G  = {}
best_codon  = {}

print("##delta G for all codon-pairs")
for i in index:
    for j in range(len(AAlist[i])):
        a = AAlist[i][j]
        c = trips[i]
        if a not in best_U:
            best_U[a]     = Utrip_random[i]
            best_G[a]     = G_random[i]
            if c[11] == a:
                best_codon[a] = c[:3]
            else:
                best_codon[a] = c[8:5:-1]
                
        else:
            if Utrip_random[i] < best_U[a]:
                best_U[a] = Utrip_random[i]
                best_G[a] = G_random[i]
                if c[11] == a:
                    best_codon[a] = c[:3]
                else:
                    best_codon[a] = c[8:5:-1]
                
#    print ("%s & %.3f & %.3f & %.3f & %.3f \\\\" %\
#       (trips[i], Utrip[i], Udelt[i], Utrip_random[i], G_random[i]))
    reparsed = trips[i][:3]+r'$\cdot$'+trips[i][8:5:-1]+" "+\
                trips[i][11:13]+trips[i][16:]
    print ("%s & %.3f & %.5f  \\\\" %\
           (reparsed, Udelt[i], np.exp(-1.*G_random[i]*beta)) )

##build a table like the one in PMC3293468
###3D list of cells, by 1 2 3 base of codon
cells = {
    'T':{ 'T':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'A':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'G':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'C':{'T':{},  'A':{}, 'G':{},  'C':{}}  },
    'C':{ 'T':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'A':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'G':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'C':{'T':{},  'A':{}, 'G':{},  'C':{}}  },
    'A':{ 'T':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'A':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'G':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'C':{'T':{},  'A':{}, 'G':{},  'C':{}}  },
    'G':{ 'T':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'A':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'G':{'T':{},  'A':{}, 'G':{},  'C':{}},\
          'C':{'T':{},  'A':{}, 'G':{},  'C':{}}  },
}

grey_fm=r'\cellcolor[gray]{%.3f} %s'
for i in index:
    cod1 = trips[i][:3]
    cod2 = trips[i][8:5:-1]
    aa1  = trips[i][11:13]
    aa2  = trips[i][16:]
    greyLevel = min(0.3+ 0.7*(np.exp(-1.*G_random[i]*beta)/0.4), 1.0)
    cells[cod1[0]][cod1[1]][cod1[2]] = grey_fm % (greyLevel, aa1)
    cells[cod2[0]][cod2[1]][cod2[2]] = grey_fm % (greyLevel, aa2)

print(r'''\begin{tabular}{  l | l | c | c | c | c |  c | l}
\multicolumn{5}{c}{Base 3}\\
\cline{3-6}
\multicolumn{2}{c|}{}&
\multicolumn{1}{c|}{T}&
\multicolumn{1}{c|}{C}&
\multicolumn{1}{c|}{A}&
\multicolumn{1}{c|}{G}\\
''')

for row in ['T', 'C', 'A', 'G']:
    print(r'\cline{2-7}')

    if row == 'T':
        print(r'''
        \multirow{16}{*}{\begin{sideways}Base 1\end{sideways}}  
               ''')
    
    for subRow in ['T', 'C', 'A', 'G']:
        if subRow == 'A':
            rowHead = '&%s' % row
        else:
            rowHead = '& '
        print("%s  & " % rowHead, end='')##python3, print supressing newline.
        for col in ['T', 'C', 'A', 'G']:
            ## middle letter is least important: so is subrow.
            cod   = row+subRow+col
            cell  = cells[cod[0]][cod[1]][cod[2]]
            print("%s &" % cell, end='') 
        
        if row == 'T' and subRow == 'T':
            print(r'''%s & \multirow{16}{*}{\begin{turn}{270}Base 2\end{turn}} 
                 \\''' % subRow, end='\n')
        else:
            print("  %s\\\\" % subRow, end='\n')
print(r'\cline{2-7}')


print("##RESULT!!!:  All phase one amino acids (except I, which synthesized from T) have an available triplet energy 0 or less.")

print("Phase1: G D E A V S T I L P")

ulist = []
alist = []
for a in best_U.keys():
    alist.append(a)
    ulist.append(best_U[a])
index = np.argsort(np.array(ulist))

for i in index:
    a = alist[i]
    if a in phase_one_aa:
        if a in phase_zero_aa:
                phase1 = "** "
        else:
                phase1 = "*  "
    else:
                phase1 = ".  "
    print("%s%s & %s & %.3f & %.4f pA: %f pB: %f pA-pB: %f\\\\" %\
 (a, phase1, best_codon[a], best_U[a], np.exp(-1.*best_G[a]*beta), prop_alpha[a], prop_beta[a], prop_alpha[a]- prop_beta[a]))

##now rank all codons
