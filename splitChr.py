#! /usr/bin/env python

import sys

with open(sys.argv[1]) as f:
  for line in f:
    if 'chrUn' not in line:
      line = line.strip()
      fields = line.split()
      chrm, length = fields[0], int(fields[1])
      if length > 5000000:
        third = length / 3
        third = int(third)
        print(f"{chrm}:0-{third}")
        print(f"{chrm}:{third+1}-{third*2}")
        print(f"{chrm}:{third*2+1}-{length}")
      else:
        print(f"{chrm}:0-{length}")
