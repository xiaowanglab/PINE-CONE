# -*- coding: utf-8 -*-
"""
@author: Kylie Standage-Beier
"""

"""To import PINE-CONE script from local directory"""
from pinecone import design_pegs #if script is renamed correct this statement

""" INPUT PARAMETERS:
organism: 'null','string' 'text','plasmid','hg38','sacCer3','mm10','danRer7','ce11','dm6'
    null = nothing 
    string, provide DNA sequence in 'file_dna' with at least 210 Bp up/downstream of edit position
    text, will retrieve sequence from txt file defined as file path in 'file_dna'
    plasmid, will retrieve txt file and treat DNA topology as circular. Path defined in 'file_dna'
    reference genomes: 'hg38', 'sacCer3','mm10','danRer7','ce11','dm6'. Will utilize API request to retrieve DNA sequences (slower)
design_alg: 'near' or 'specific' determine protospacer selection for peg and PE3
edit_strat: 'PE2', 'PE3' or 'PE3B'
name: 'pegRNA name'
chomosome: 'string number or other'
postion: 'string number'
edit: 'string sequence'
rt: 'string number'
pbs: 'string number'
notes: 'string notes'
file_dna: 'string of file path or DNA seq'
"""

""" Example Data in """
organism = 'hg38'
design_alg = 'specific'
edit_strat = 'PE3'
name = 'APOE24'
chromosome = '19'
position = '44908756'
edit = 'T'
rt = '24'
pbs = '13'
notes= 'APOE(R136S) A'
file_dna = ''

""" Data is entered as list of dictionaries, [{pegRNA-1},{pegRNA-2},{pegRNA-3} etc] """

data_in =  [{'organism':organism,'design_alg':design_alg,'edit_strat':edit_strat,'name':name,'chromosome':chromosome,'position':position,'edit':edit,'rt':rt,'pbs':pbs,'notes':notes,'file_dna':file_dna}]

""" Design pegRNAs is called via 'design_pegs' function. Data out is returned as list of dictionaries """
data_out = design_pegs(data_in)
print(data_out)

""" OUTPUT 'KEY:VALUE' PAIRS:
'index':assigned guide index, 
'name':pegRNA_name
'chromosome':chromosome
'position':edit position
'product':'X-to-Y' editing product
'strand':'positive' or 'negative' strand
'pegrna':pegRNA sequence
'proto_peg':pegRNA protospacer
'rt':RT sequence
'rt_length':RT length
'pbs':PBS sequence
'pbs_length':PBS length
'proto_pe3':PE3 or 3B protospacer
'pe3_distance':Cut distance (peg to PE3)
'peg_guide_top':Oligo pegRNA guide top
'peg_guide_bot':Oligo pegRNA guide bottom
'peg_ext_top':Oligo pegRNA extenstion top
'peg_ext_bot':Oligo pegRNA extenstion bottom
'pe3_top':PE3 or 3B guide top
'pe3_bot':PE3 or 3B guide bottom
'primer_f':estimated PCR primer forward
'primer_r':estimated PCR primer reverse
'notes':notes appended poly-U track warnings 
'wt_seq':WT Edit position +/- 30 Bp
'edit_seq':Edited DNA +/- 30 Bp
"""