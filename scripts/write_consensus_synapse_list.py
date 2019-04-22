import sys
sys.path.append('./volumetric_analysis')
sys.path.append('.')
from lxml import etree
from tqdm import tqdm

from mat_loader import MatLoader
import db
import aux

def scrub_left_chemical_synapses(cur,C,cells,start=0,end=None):
    data = {}

    for cname in tqdm(cells,desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_presynapse_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            post,neigh,sect = [],[],[]
            for s in syn:
                sect.append(s[0])
                adj = db.mine.get_object_adjacency(cur,s[1])
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(cname,n): continue
                    if n in neigh: continue
                    neigh.append(n)
 
                for p in s[2]:
                    pname = db.mine.get_object_contin_name(cur,p)
                    if not pname: continue
                    if not C.C.has_edge(cname,pname): continue 
                    if pname in post: continue
                    if pname not in neigh: continue
                    post.append(pname)

            if not post: continue
            data[cname][cont] = {'partners':sorted(post),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data

def scrub_right_chemical_synapses(cur,C,cells,lrmap,start=0,end=None):
    data = {}
    for cname in tqdm(cells,desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_presynapse_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            post,neigh,sect = [],[],[]
            for s in syn:
                sect.append(s[0])
                adj = db.mine.get_object_adjacency(cur,s[1])
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(lrmap[cname],lrmap[n]): continue
                    if n in neigh: continue
                    neigh.append(n)
 
                for p in s[2]:
                    pname = db.mine.get_object_contin_name(cur,p)
                    if not pname: continue
                    if pname not in lrmap: continue
                    if not C.C.has_edge(lrmap[cname],lrmap[pname]): continue 
                    if pname in post: continue
                    if pname not in neigh: continue
                    post.append(pname)

            if not post: continue
            data[cname][cont] = {'partners':sorted(post),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data

def scrub_left_chemical_post_synapses(cur,C,cells,start=0,end=None):
    data = {}

    for cname in tqdm(cells,desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_postsynapse_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            pre,neigh,sect = [],[],[]
            for s in syn:
                sect.append(s[0])
                pobj = None
                for p in s[2]:
                    pname = db.mine.get_object_contin_name(cur,p)
                    if not pname: continue
                    if pname == cname: 
                        pobj = p
                        break
                if not pobj: continue
                adj = db.mine.get_object_adjacency(cur,pobj)
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(cname,n): continue
                    if n in neigh: continue
                    neigh.append(n)
 
                
                pname = db.mine.get_object_contin_name(cur,s[1])
                if not pname: continue
                if pname in pre: continue
                if not C.C.has_edge(pname,cname): continue
                if pname == cname: continue
                if pname not in neigh: continue
                pre.append(pname)

            if not pre: continue
            data[cname][cont] = {'partners':sorted(pre),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data

def scrub_right_chemical_post_synapses(cur,C,cells,lrmap,start=0,end=None):
    data = {}

    for cname in tqdm(cells,desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_postsynapse_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            pre,neigh,sect = [],[],[]
            for s in syn:
                sect.append(s[0])
                pobj = None
                for p in s[2]:
                    pname = db.mine.get_object_contin_name(cur,p)
                    if not pname: continue
                    if pname == cname: 
                        pobj = p
                        break
                if not pobj: continue
                adj = db.mine.get_object_adjacency(cur,pobj)
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(lrmap[cname],n): continue
                    if n in neigh: continue
                    neigh.append(n)
                
                pname = db.mine.get_object_contin_name(cur,s[1])
                if not pname: continue
                if pname in pre: continue
                if pname not in lrmap: continue
                if not C.C.has_edge(lrmap[pname],lrmap[cname]): continue
                if pname == cname: continue
                if pname not in neigh: continue
                pre.append(pname)

            if not pre: continue
            data[cname][cont] = {'partners':sorted(pre),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data


def scrub_left_gapjunction_synapses(cur,C,cells,start=0,end=None):
    data = {}
    for cname in tqdm(cells, desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_gapjunction_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            post,neigh,sect = [],[],[]
            for s in syn:
                preobj,postobj = s[1],s[2][0]
                if not db.mine.get_object_contin_name(cur,s[1]) == cname:
                    preobj,postobj = postobj,preobj
                
                adj = db.mine.get_object_adjacency(cur,preobj)
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(cname,n): continue
                    if n in neigh: continue
                    neigh.append(n)
                
                pname = db.mine.get_object_contin_name(cur,postobj)
                if not pname: continue
                if pname in post: continue
                if not C.E.has_edge(cname,pname): continue 
                if pname in post: continue
                if pname not in neigh: continue
                post.append(pname)

            if not post: continue
            data[cname][cont] = {'partners':sorted(post),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data

def scrub_right_gapjunction_synapses(cur,C,cells,lrmap,start=0,end=None):
    data = {}
    for cname in tqdm(cells,desc="Cells processed"):
        data[cname] = {}
        contins = db.mine.get_gapjunction_contins(cur,cname,start=start,end=end)
        for cont in contins:
            syn = db.mine.get_synapse_data_by_contin(cur,cont)
            if not syn: continue
            post,neigh,sect = [],[],[]
            for s in syn:
                preobj,postobj = s[1],s[2][0]
                if not db.mine.get_object_contin_name(cur,s[1]) == cname:
                    preobj,postobj = postobj,preobj
                
                adj = db.mine.get_object_adjacency(cur,preobj)
                if not adj: continue
                for n in adj:
                    if not C.A.has_edge(lrmap[cname],lrmap[n]): continue
                    if n in neigh: continue
                    neigh.append(n)
                
                pname = db.mine.get_object_contin_name(cur,postobj)
                if not pname: continue
                if pname in post: continue
                if pname not in lrmap: continue
                if not C.E.has_edge(lrmap[cname],lrmap[pname]): continue 
                if pname in post: continue
                if pname not in neigh: continue
                post.append(pname)

            if not post: continue
            data[cname][cont] = {'partners':sorted(post),
                                'neighbors':sorted(neigh),
                                'num_sections':len(syn),
                                'sections':sect}
    return data



def convert_to_xml(data):
    root = etree.Element('conserved_synapses')
    for cell in data:
        xcell = etree.SubElement(root,'cell')
        xcell.set('name',cell)
        for cont in data[cell]:
            xcont = etree.SubElement(xcell,'contin')
            xcont.set('id',cont)
            xnum = etree.SubElement(xcont,'num_sections')
            xnum.text = str(data[cell][cont]['num_sections'])
            xsect = etree.SubElement(xcont,'sections')
            for s in data[cell][cont]['sections']:
                xname = etree.SubElement(xsect,'name')
                xname.text = s
            xpartner = etree.SubElement(xcont,'partners')
            for p in data[cell][cont]['partners']:
                xname = etree.SubElement(xpartner,'name')
                xname.text = p
            xneighbor = etree.SubElement(xcont,'neighbors')
            for n in data[cell][cont]['neighbors']:
                xname = etree.SubElement(xneighbor,'name')
                xname.text = n
            xcell.append(xcont)
        root.append(xcell)

    tree = etree.ElementTree(root)
    return tree      

DEG = [3,4]
DB = ['N2U','JSH']
PROCESS = [0,0,1,1,0,0]

if __name__=="__main__":

    for deg in DEG:
        M = MatLoader()
        M.load_left()
        M.load_right()
        M.load_lrmap()
        C = M.load_consensus_graphs(deg)

        for _db in DB:
            start,end = 0,325
            if _db == 'JSH':
                start,end = 0,425

            clout = './mat/consensus_synapses/%s_chem_synapses_left_deg%d.xml' %(_db,deg)
            crout = './mat/consensus_synapses/%s_chem_synapses_right_deg%d.xml' %(_db,deg)
            cplout = './mat/consensus_synapses/%s_chem_post_synapses_left_deg%d.xml' %(_db,deg)
            cprout = './mat/consensus_synapses/%s_chem_post_synapses_right_deg%d.xml' %(_db,deg)
            elout = './mat/consensus_synapses/%s_elec_synapses_left_deg%d.xml' %(_db,deg)            
            erout = './mat/consensus_synapses/%s_elec_synapses_right_deg%d.xml' %(_db,deg)

            con = db.connect.default(_db)
            cur = con.cursor()
            
            if PROCESS[0]:
                print('Left chemical, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_left_chemical_synapses(cur,C,M.left,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(clout,'wb') as fout: fout.write(xml_out)
            
            if PROCESS[1]:
                print('Right chemical, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_right_chemical_synapses(cur,C,M.right,M.lrmap,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(crout,'wb') as fout: fout.write(xml_out)
            
            if PROCESS[2]:
                print('Left chemical post, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_left_chemical_post_synapses(cur,C,M.left,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(cplout,'wb') as fout: fout.write(xml_out)
            
            if PROCESS[3]:
                print('Right chemical post, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_right_chemical_post_synapses(cur,C,M.right,M.lrmap,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(cprout,'wb') as fout: fout.write(xml_out)

            if PROCESS[4]:
                print('Right gap junction, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_left_gapjunction_synapses(cur,C,M.left,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(elout,'wb') as fout: fout.write(xml_out)

            if PROCESS[5]:
                print('Right gap junction, DB: %s, DEG: %d' %(_db,deg))
                data = scrub_right_gapjunction_synapses(cur,C,M.right,M.lrmap,start=start,end=end)
                tree = convert_to_xml(data)
                xml_out = etree.tostring(tree,pretty_print=False)
                with open(erout,'wb') as fout: fout.write(xml_out)

            con. close()
