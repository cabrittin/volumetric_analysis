import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')

from mat_loader import MatLoader
import db
import aux

def scrub_synapses(G,synapses):
    scrubbed = []
    for [pre,_post,sect,contin,series] in synapses:
        if not C.C.has_node(pre): continue
        _post = _post.split(',')
        post = [p for p in _post if C.C.has_edge(pre,p)]
        tmp = [pre,post,sect,contin,series]
        if not post: continue
        scrubbed.append(tmp)
    return scrubbed


if __name__=="__main__":
    
    for deg in  [3,4]:
        M = MatLoader()
        C = M.load_consensus_graphs(deg)

        for _db in  ['N2U','JSH']:
            start,end = 0,325
            if _db == 'JSH':
                start,end = 0,425

            cout = './mat/consensus_synapses/%s_chem_synapses_deg%d.csv' %(_db,deg)
            eout = './mat/consensus_synapses/%s_elec_synapses_deg%d.csv' %(_db,deg)
    
            con = db.connect.default(_db)
            cur = con.cursor()

            synapses = db.mine.get_synapse_data(cur,'chemical',start=start,end=end)
            scrub = scrub_synapses(C.C,synapses)
            aux.write.from_list(cout,scrub)
    
            synapses = db.mine.get_synapse_data(cur,'electrical',start=start,end=end)
            scrub = scrub_synapses(C.E,synapses)
            aux.write.from_list(eout,scrub)
    
            con. close()
