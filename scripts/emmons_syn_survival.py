import sys
sys.path.append('./volumetric_analysis')
import matplotlib.pyplot as plt
import aux

HERM = "#EF24A4"
MALE = "#A845EF"

herm = aux.read.into_list2('./data/herm_synapses.csv',delimiter='\t')
male = aux.read.into_list2('./data/male_synapses.csv',delimiter='\t')

hc,he = [],[]
for h in herm[1:]:
    if h[0] not in ['N2U','JSE']: continue
    c = 1
    if h[2] == 'Not NR': c = 2  
    if h[5] == 'chemical': hc.append(int(h[6])*c)
    if h[5] == 'electrical': he.append(int(h[6])*c)

mc,me = [],[]
for m in male[1:]:
    if m[5] == 'chemical': mc.append(int(m[6]))
    if m[5] == 'electrical' : me.append(int(m[6]))

    

plt.figure(figsize=(10,12))
plt.hist(hc,bins=100,range=(1,100),histtype='step',cumulative=-1,density=True,color=HERM,linewidth=2,label="Hermaphrodite (n=%d)"%len(hc))
plt.hist(mc,bins=100,range=(1,100),histtype='step',cumulative=-1,density=True,color=MALE,linewidth=2,label="Male (n=%d)"%len(mc))
plt.yscale('log')
plt.xlim([1,55])
plt.title('Chemical synapses',fontsize=34)
plt.ylabel('Survival function',fontsize=28)
plt.xlabel('Synapse size (# EM sections)',fontsize=28)
plt.legend(fontsize=18)
plt.savefig('chemical_synapses.png')

plt.figure(figsize=(10,12))
plt.hist(he,bins=100,range=(1,100),histtype='step',cumulative=-1,density=True,color=HERM,linewidth=2,label="Hermaphrodite (n=%d)"%len(he))
plt.hist(me,bins=100,range=(1,100),histtype='step',cumulative=-1,density=True,color=MALE,linewidth=2,label="Male (n=%d)"%len(me))
plt.yscale('log')
plt.xlim([1,55])
plt.title('Gap junctions',fontsize=34)
plt.ylabel('Survival function',fontsize=28)
plt.xlabel('Synapse size (# EM sections)',fontsize=28)
plt.legend(fontsize=18)
plt.savefig('gap_junctions.png')


plt.show()


