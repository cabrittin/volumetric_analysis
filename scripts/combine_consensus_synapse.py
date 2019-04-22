"""
combine_consensus_synapse.py

Combines consensus synapse data

"""

import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
from lxml import etree

from connectome import consensus
from mat_loader import MatLoader

def combine_data(data0,data1):
    for cell in data1:
        if cell in data0:
            for cont in data1[cell]:
                if cont in data0[cell]:
                    newcont = cont
                    while newcont in data0[cell]:
                        newcont = '_' + newcont
                    data0[cell][newcont] = data1[cell][cont]
                else:
                    data0[cell][cont] = data1[cell][cont]
        else:
            data0[cell] = data1[cell]

    return data0

def map_right_left(rdata,lrmap):
    ldata = {}
    for rcell in rdata:
        lcell = lrmap[rcell]
        ldata[lcell] = {}
        for cont in rdata[rcell]:
            lpartners = [lrmap[rp] for rp in rdata[rcell][cont]['partners']]
            lneighbors = [lrmap[rn] for rn in rdata[rcell][cont]['neighbors']]
            ldata[lcell][cont] = {'num_sections': rdata[rcell][cont]['num_sections'],
                                'sections' : rdata[rcell][cont]['sections'],
                                'partners' : lpartners,
                                'neighbors' : lneighbors}
    return ldata

din = './mat/consensus_synapses/'
chem_left = ['N2U_chem_synapses_left_deg%d.xml',
            'JSH_chem_synapses_left_deg%d.xml']
chem_right = ['N2U_chem_synapses_right_deg%d.xml',
            'JSH_chem_synapses_right_deg%d.xml']
post_left = ['N2U_chem_post_synapses_left_deg%d.xml',
            'JSH_chem_post_synapses_left_deg%d.xml']
post_right = ['N2U_chem_post_synapses_right_deg%d.xml',
            'JSH_chem_post_synapses_right_deg%d.xml']
gap_left = ['N2U_elec_synapses_left_deg%d.xml',
            'JSH_elec_synapses_left_deg%d.xml']
gap_right = ['N2U_elec_synapses_right_deg%d.xml',
            'JSH_elec_synapses_right_deg%d.xml']

cout =  din + 'all_chem_deg%d.xml'
pout = din + 'all_chem_post_deg%d.xml'
eout = din + 'all_elec_deg%d.xml'

DEG = [1,2,3,4]
if __name__=="__main__":
    M = MatLoader()
    M.load_lrmap()
    for deg in DEG:
        data = consensus.convert_xml_to_synapse(din + chem_left[0]%deg)
        tmp = consensus.convert_xml_to_synapse(din + chem_left[1]%deg)
        data = combine_data(data,tmp)

        tmp = consensus.convert_xml_to_synapse(din + chem_right[0]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)
        
        tmp = consensus.convert_xml_to_synapse(din + chem_right[1]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)

        tree = consensus.convert_synapse_to_xml(data)
        xml_out = etree.tostring(tree,pretty_print=False)
        with open(cout%deg,'wb') as fout: fout.write(xml_out)

        data = consensus.convert_xml_to_synapse(din + post_left[0]%deg)
        tmp = consensus.convert_xml_to_synapse(din + post_left[1]%deg)
        data = combine_data(data,tmp)

        tmp = consensus.convert_xml_to_synapse(din + post_right[0]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)

        tmp = consensus.convert_xml_to_synapse(din + post_right[1]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)

        tree = consensus.convert_synapse_to_xml(data)
        xml_out = etree.tostring(tree,pretty_print=False)
        with open(pout%deg,'wb') as fout: fout.write(xml_out)

        data = consensus.convert_xml_to_synapse(din + gap_left[0]%deg)
        tmp = consensus.convert_xml_to_synapse(din + gap_left[1]%deg)
        data = combine_data(data,tmp)

        tmp = consensus.convert_xml_to_synapse(din + gap_right[0]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)
        
        tmp = consensus.convert_xml_to_synapse(din + gap_right[1]%deg)
        tmp = map_right_left(tmp,M.lrmap)
        data = combine_data(data,tmp)

        tree = consensus.convert_synapse_to_xml(data)
        xml_out = etree.tostring(tree,pretty_print=False)
        with open(eout%deg,'wb') as fout: fout.write(xml_out)



         
