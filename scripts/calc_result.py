#!/usr/bin/env python3
import sys
fn = sys.argv[1]
with open(fn, 'r') as f:
  read2id = dict()
  for line in f:
    if line.startswith('@taxid'):
      read, _, taxid = line.split(' ')
      taxid = taxid.split('/')[-1].split('.')[0]
      #genus_id = get_genusid(taxid.split('_')[-1].split('.')[0]) #
      read = read.split('/')[0]
      #if (read in read2id) and (taxid != read2id[read][0]):
      #  print('error:')
      if read in read2id:
        if taxid != read2id[read]:
          #print('find paired', read, taxid)
          read2id[read] = -1
      else:
        read2id[read] = taxid

count = 0
for k, v in read2id.items():
  if v != -1:
    count += 1
    print(k, v)

print(len(read2id), count)
