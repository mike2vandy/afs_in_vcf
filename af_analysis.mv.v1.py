#! /usr/bin/env python

import sys, vcf
from collections import Counter

def reduceCSQ(csqs):
  #columns: variant,consequence,impact,gene,transcriptID,biotype,exon,hgvsc,hgvsp,cDNA_pos,cds_pos,aa_pos,aas,codons,variant_class,proteinID
  mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 27, 28, 30]
  #mainCols = [0, 1, 2, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 21, 30]
  updated = []
  for record in csqs:
    record = record.split('|')
    result = ','.join([record[i] for i in mainCols])
    updated.append(result)
  return(updated)

def calcAF(groupList):
  #calculate MAF
  alleles = [allele for genotype in groupList for allele in genotype.split('/') if allele != '.']
  altAF = 1 - (alleles.count('0') / len(alleles))
  #maf = round(min(altAF, 1 - altAF), 6)
  
  #Count Genotypes
  genotypes = Counter(groupList)
  homRef = genotypes.get('0/0', 0)
  miss = genotypes.get('./.', 0)
  het = 0
  homVar = 0
  for genos, count in genotypes.items():
    if genos != '0/0' and genos != './.':
      geno = genos.split('/')
      if geno[0] == geno[1]:
        homVar += count
      else:
        het += count
  return([homRef, het, homVar, miss, altAF]) 

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
  AF = line.INFO.get('AF')
  AF = sum(AF) 
  if AF < 0.5:
    major, minor, maf = ref, alt, AF
  else:
    major, minor, maf = alt, ref, 1-AF
  maf = round(maf, 6)
  first = ','.join(map(str, [chrm, pos, ref, alt, major, minor, maf]))
  ctl, risk, aff, main = [], [], [], []
  for i in line.samples:
    sample, gt = i.sample, i['GT'].replace('|', '/')
    main.append(gt)
    group = dogs.get(sample, None)
    if group == 'ctl':
      ctl.append(gt)
    elif group == 'risk':
      risk.append(gt)
    elif group == 'aff':
      aff.append(gt)
 
  countsRisk = ','.join(map(str, calcAF(risk)))
  countsAff = ','.join(map(str, calcAF(aff)))
  countsCtl = ','.join(map(str, calcAF(ctl)))
  countsAll = ','.join(map(str,calcAF(main)[:4]))
  csq = reduceCSQ(line.INFO.get('CSQ', []))
  for consequence in csq:
    print(f"{first},{countsAff},{countsRisk},{countsCtl},{countsAll},{consequence}")
  
     
