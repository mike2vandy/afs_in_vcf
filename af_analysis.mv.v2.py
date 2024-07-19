#!/usr/bin/env python

import os
import sys
from collections import Counter

import vcf

def reduceCSQ(csqs):
    #columns: variant,consequence,impact,gene,transcriptID,biotype,exon,hgvsc,hgvsp,cDNA_pos,cds_pos,aa_pos,aas,codons,variant_class,proteinID
    #mainCols = [0, 1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 21, 30]
    mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 27, 28, 30]
    updated = []
    for record in csqs:
        record = record.split('|')
        if '&' in record[-1]:
           #print('HERE IS ONE')
            record[-1] = "&".join(sorted(set(record[-1].split('&')), key=record[-1].split('&').index))
        result = ','.join([record[i] for i in mainCols])
        updated.append(result)
    
    return(updated)

def majority_element(nums):
    return '1' if sum(int(num) for num in nums) > len(nums) / 2 else '0'

def calcAF(groupList, minorAlt=True):
    #calculate MAF
    alleles = [allele for genotype in groupList for allele in genotype.split('/') if allele != '.']
    alt_count = '1' 
    if not minorAlt:
        alt_count = '0'
    try:
       #altAF = 1 - (alleles.count(alt_count) / len(alleles))
       #maf = round(min(altAF, 1 - altAF), 6)
        maf = round(alleles.count(alt_count)/len(alleles), 6)
    except ZeroDivisionError:
        maf = 'NA'
  
    #Count Genotypes
    genotypes = Counter(groupList)
    homRef = genotypes.get('0/0', 0)
    miss = genotypes.get('./.', genotypes.get('.',0))
    het = 0
    homVar = 0
    for genos, count in genotypes.items():
        if genos != '0/0' and genos != './.' and genos != '.':
            geno = genos.split('/')
            if geno[0] == geno[1]:
                homVar += count
            else:
                het += count
    
    return([homRef, het, homVar, miss, maf]) 

#parse dog list  
dogs = {}
for i in sys.argv[2:]:
    group = i.split('.')[0]
    with open(i) as f: 
        for sample in f:
            sample = sample.strip()
            dogs[sample] = group

#read vcf
vcfR = vcf.Reader(filename = sys.argv[1])

for line in vcfR:
    chrm, pos, ref, alt = line.CHROM, line.POS, line.REF, '/'.join(map(str, line.ALT))
    if line.FILTER:
      filt = line.FILTER[0]
    else:
      filt = 'PASS'
    AF = line.INFO.get('AF')
    AF = sum(AF)
    minor_alt = True
    if AF < 0.5:
        major, minor, maf = ref, alt, AF
    else:
        major, minor, maf = alt, ref, 1-AF
        minor_alt = False
    maf = round(maf, 6)
    first = ','.join(map(str, [chrm, pos, filt, ref, alt, major, minor, maf]))
   #ctl, risk, aff, main = [], [], [], []
    main, tvd, lab, sick, ctrl = [], [], [], [], []
    for i in line.samples:
        sample, gt = i.sample, i['GT'].replace('|', '/') 
        main.append(gt)
        group = os.path.basename(dogs.get(sample, "")) if dogs.get(sample) is not None else None 

        if group == 'ctlDogs':
            ctrl.append(gt)
        elif group == 'tvds':
            tvd.append(gt)
        elif group == 'sickHearts':
            sick.append(gt)
        elif group == 'otherLabs':
            lab.append(gt)


    countsTvd = ','.join(map(str, calcAF(tvd, minor_alt)))
    countsLab = ','.join(map(str, calcAF(lab, minor_alt)))
    countsSick = ','.join(map(str, calcAF(sick, minor_alt)))
    countsCtl = ','.join(map(str, calcAF(ctrl, minor_alt)))
    countsAll = ','.join(map(str,calcAF(main, minor_alt)[:4]))
    csq = reduceCSQ(line.INFO.get('CSQ', []))
    for consequence in csq:
        print(f"{first},{countsTvd},{countsLab},{countsSick},{countsCtl},{countsAll},{consequence}")
  
     
