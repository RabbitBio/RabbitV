
target = [
  '10804',
  '1089136',
  '11053',
  '1136133',
  '11801',
  '1327964',
  '1544901',
  '198503',
  '37138',
  '915425',
]
#target="""
#1036673
#1042876
#1053692
#1104326
#1117943
#1123863
#1133568
#1159202
#1184253
#1218356
#122586
#1328311
#226185
#243276
#262728
#290339
#299768
#300852
#347495
#354242
#367737
#373384
#395492
#402882
#418127
#434271
#440085
#441952
#449216
#536232
#552536
#568706
#573059
#573236
#768492
#862962
#889738
#930943
#990315
#998820"""
#target = target.split('\n')

import glob
import os
import subprocess
from multiprocessing import Process
import get_tax_uplvl

def grep_from_dir(name, file_dir):
  for filename in glob.glob(file_dir + '*.fna'):
    with open(filename, 'r') as f:
      for line in f:
        if line.find(name) != -1:
          #print('find:', line)
          return filename
        
def process_func(name, taxid):
  file_path = grep_from_dir(name, '/home/ssd/viral/')
  return file_path
  #cmd = ['cp', file_path, os.path.join('./strain_excluded', taxid+'.fa')]
  #subprocess.check_output(cmd)

nodes_p_r_dict = get_tax_uplvl.nodes2dict(get_tax_uplvl.nodes_file)
with open("./names.dmp") as f:
  ps = []
  for line in f:
    taxid, name, _, comm, _ = [x.strip() for x in line.split('|')]
    if taxid in target and comm == 'scientific name':
      file_path = process_func(name, taxid)
      genus_taxid = get_tax_uplvl.find_parent_in_lvl(taxid, nodes_p_r_dict=nodes_p_r_dict, taxonomy_lvl = 'genus')
      #print("writing file: {} to file: {}".format(file_path, genus_taxid))
      if not genus_taxid:
        print('error to find the {} level data of {}, we use the taxid as filename!'.format(taxonomy_lvl, taxid))
        print('line: ', line)
        exit(-1)
        filename = taxid + '.fa'
      else:
        filename = genus_taxid + '.fa'

      #print(taxid, name, ': cat', file_path, '>>', genus_taxid+'.fa')
      print(taxid, name, ' cp', file_path, taxid+'.fa')
      #ref_file = open(os.path.join('./strain_excluded_genus/', filename), 'w+')
      #ref_file.write()
      #ref_file.close()


  print('done')
