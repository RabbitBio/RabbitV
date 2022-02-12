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
    #print("find: ", input_taxid)
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


def taxid2parent(input_taxid_list, taxonomy_lvl='genus', nodes_file='./nodes.dmp'):
    return_list = []
    nodes_p_r_dict = nodes2dict(nodes_file)
    # modified input list
    if type(input_taxid_list) is list:
        input_taxid_list = input_taxid_list
    else:
        input_taxid_list = eval(input_taxid_list)
    # find taxid parent 
    for input_taxid in input_taxid_list:
        if input_taxid not in nodes_p_r_dict.keys():
            print(input_taxid, '\ttaxid_error')
            return_list.append([input_taxid, 'taxid_error'])
        else:
            parent_taxid = find_parent_in_lvl(str(input_taxid), nodes_p_r_dict=nodes_p_r_dict, taxonomy_lvl='genus')
            print(f'{input_taxid}\t{parent_taxid}')
            return_list.append([input_taxid, parent_taxid])
    return return_list


if __name__ == '__main__':
  #target = [
  #  '10804',
  #  '1089136',
  #  '11053',
  #  '1136133',
  #  '11801',
  #  '1327964',
  #  '1544901',
  #  '198503',
  #  '37138',
  #  '915425',
  #]
  target="""
  1036673
  1042876
  1053692
  1104326
  1117943
  1123863
  1133568
  1159202
  1184253
  1218356
  122586
  1328311
  226185
  243276
  262728
  290339
  299768
  300852
  347495
  354242
  367737
  373384
  395492
  402882
  418127
  434271
  440085
  441952
  449216
  536232
  552536
  568706
  573059
  573236
  768492
  862962
  889738
  930943
  990315
  998820"""
  target = [x.strip() for x in target.split('\n')]
  input_taxid_list = target  #[9606,35]
  parent_taxid_list = taxid2parent(input_taxid_list)
