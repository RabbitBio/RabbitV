#!/usr/bin/env python3
import os

import sys

nodes_file = './nodes.dmp'


def nodes2dict(nodes_file):
  nodes_p_r_dict = {}
  node_data_list = []
  with open(nodes_file) as f_nodes:
    #node_data_list = [i.split('|') for i in f_nodes.read().strip().split('\n')]
    for line in f_nodes:
      node_data_list.append([x.strip() for x in line.split('|')])

  for line_list in node_data_list:
    taxid = line_list[0]
    parent_taxid = line_list[1]
    rank = line_list[2]
    nodes_p_r_dict[taxid] = {'parent_taxid': parent_taxid,
                              'rank': rank}
  return nodes_p_r_dict


def find_parent_in_lvl(input_taxid, nodes_p_r_dict, taxonomy_lvl='genus', up_count=0):
    up_count += 1
    if up_count >= 20:
        return False
    parent_taxid = nodes_p_r_dict[input_taxid]['parent_taxid']
    rank = nodes_p_r_dict[input_taxid]['rank']
    # print(rank, taxonomy_lvl)
    if rank == taxonomy_lvl:
        #print('find:', input_taxid, rank)
        return input_taxid
    else:
        return find_parent_in_lvl(parent_taxid, nodes_p_r_dict, taxonomy_lvl, up_count)

def mk_ref_dif(fpath, taxonomy_lvl):
  dir_name = os.path.join('.', fpath.split('/')[-1].split('.')[0]+'_'+taxonomy_lvl) 
  if not os.path.exists(dir_name):
    os.mkdir(dir_name)
  else:
    print("error: dir exist!")
    exit(-1)
  return dir_name

def split_by_species():
  reference_all = '/home/user_home/public/kraken2-reference-data/strain_excluded.fna'
  taxonomy_lvl = 'species'
  dir_name = mk_ref_dif(reference_all, taxonomy_lvl) #make the dir by taxnomy level
  with open(reference_all, 'r') as f:
    ref_file = None
    for line in f:
      if line.startswith('>'):
        #>gi|170018061|ref|NC_010468.1|kraken:taxid|481805| Escherichia coli ATCC 8739 chromosome, complete genome
        taxid = line.split('|')[5]
        filename = os.path.join(dir_name, taxid + '.fa')

        if ref_file is None :
          ref_file = open(filename, 'a+')
        if filename != ref_file.name: #open a new file
          ref_file.close()
          ref_file = open(filename, 'a+')
          ref_file.write(line)
      else:
        ref_file.write(line)
    ref_file.close()

def split_by_genus():
  reference_all = '/home/user_home/public/kraken2-reference-data/strain_excluded.fna'
  taxonomy_lvl='genus'
  
  dir_name = mk_ref_dif(reference_all, taxonomy_lvl) #make the dir by taxnomy level
  nodes_p_r_dict = nodes2dict(nodes_file)  # make 'nodeid: {parentid, rank}' dict

  with open(reference_all, 'r') as f:
    ref_file = None
    for line in f:
      if line.startswith('>'):
        #>gi|170018061|ref|NC_010468.1|kraken:taxid|481805| Escherichia coli ATCC 8739 chromosome, complete genome
        taxid = line.split('|')[5]
        genus_taxid = find_parent_in_lvl(taxid, nodes_p_r_dict=nodes_p_r_dict, taxonomy_lvl=taxonomy_lvl)
        filename = ''
        if not genus_taxid:
          print('error to find the {} level data of {}, we use the taxid as filename!'.format(taxonomy_lvl, taxid))
          print('line: ', line)
          filename = taxid + '.fa'
        else:
          filename = genus_taxid + '.fa'

        if ref_file is None :
          ref_file = open(os.path.join(dir_name, filename), 'a+')
        if filename != ref_file.name: #open a new file
          ref_file.close()
          ref_file = open(os.path.join(dir_name, filename), 'a+')
          ref_file.write(line)
      else:
        ref_file.write(line)
    ref_file.close()

if __name__ == '__main__':
  split_by_species()
