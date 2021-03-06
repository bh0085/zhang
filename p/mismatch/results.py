cluc = '''75862100		93590300		93272700		48169000		50309500		54096800	
401106000		303828000		305798000		26866100		28209200		30129900	
16983000		11725400		9053830		13685800		14116800		16643800	
38474200		29856200		27934000		315057000		311363000		317580000	
92350400		100007000		90783200		52040900		48166900		56552600	
70587900		55716200		49522200		85831900		82354100		95085500	
72628500		75919700		60881300		23941400		22015600		24504000	
26164900		28161200		22801200		10355700		9491790		9444750	'''

gluc = '''263439		395252		427727		541115		593935		591172
696121		647024		672936		662691		705704		660862
719609		657810		647289		687993		644459		669953
691331		669607		643062		231829		218895		205300
185286		206475		197801		186040		172390		174610
196630		183765		171655		174240		160860		160509
157919		170391		158607		157033		150888		133329
2968120		4138900		3493480		3824990		3622580		3084780 '''

from numpy import *
import cb.utils.plots as myplots

def get_data():
    clrows = array([[ float(c) for c in l.split('\t\t')] for l in cluc.splitlines()])
    glrows = array([[ float(c) for c in l.split('\t\t')] for l in gluc.splitlines()])
    
    return {'cluc':clrows, 'gluc':glrows}


def analyze():
    f = myplots.fignum(1)
    
    gl = get_data()['gluc']
    cl = get_data()['cluc']

    ax = f.add_subplot(111)
    ax.imshow(cl / gl, aspect = 'auto', interpolation = 'nearest')
    ax.set_title('Enrichment of Cluc over control Gluc')
    path = myplots.figpath('corr_matrix.pdf')
    f.savefig(path)
    
    f.clear()
    ax = f.add_subplot(121)
    glf = gl.flatten()[:-6]
    clf = cl.flatten()[:-6]
    
    n_mm = array([ [e,e,e] for e in [0,2,2,2,3,3,3,0,2,2,2,3,3,3]], float).flatten()
    
    
    ax.set_title('cluc enrichment vs mm count')
    ax.set_xlabel('mismatch count')
    ax.set_ylabel('fold enrichment ocluc')
    ax.scatter(n_mm,clf/glf)
    pf1 = polyfit(n_mm, clf/glf, 1)
    pf2 = polyfit(n_mm, clf/glf, 2)
    
    ax.plot(polyval(pf1,[0,2,3]))
   
    path = myplots.figpath('enrichment_vs_mm.pdf')
    f.savefig(path)

    ax2 = f.add_subplot(121)
