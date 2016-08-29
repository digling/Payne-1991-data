from lingpy import *
from lingpy.algorithm import misc

# open the file and read it line by line
D = {}
def _modify(char):
    
    char = char.strip()

    #if char == '[':
    #    return '('
    #elif char == ']':
    #    return ')'
    if not char:
        return '-'
    return char

tnames = dict([(a,(b,c)) for a,b,c in csv2list('langs.csv')])

idx = 1
cogid = 1

R = {
        'PLA' : 'PAL',
        'TAV' : 'YAV',
        'JUC' : 'YUC',
        }
with open('Payne_Maipuran_alinhamentos.xlsx - Sheet1.tsv') as f:
    
    old_concept = ''
    f.readline()
    for line in f:
        print(line)
        
        cells = line.split('\t')
        note = cells[0].strip()
        concept = cells[1].strip()

        if not concept:
            raise ValueError("no concept")

        taxon = cells[2].strip().upper()
        otaxon = taxon

        if taxon == 'PROTO-MARIUPE':
            taxon = "PROTO-MAIPURE"
        
        taxon = R.get(taxon, taxon)

        if otaxon == taxon:
            otaxon = ''

        print(taxon)
        taxon_name = tnames[taxon][0]
        iso = tnames[taxon][1]
        alm = [a.strip() for a in cells[3:]]

        if old_concept != concept:
            old_concept = concept
            cogid += 1

        # we parse the alignments in a simple manner first
        ipa = ''.join([a for a in alm])
        tokens = ' '.join([a for a in alm if a and a not in '[]'])
        alignment = ' '.join([_modify(a) for a in alm])

        # original alignment
        oralm = ' '.join([a if a else '-' for a in alm])
        
        D[idx] = [taxon, taxon_name, otaxon, iso, concept, ipa, tokens, alignment,
                oralm, cogid, note]
        idx += 1
D[0] = ['taxon', 'taxon_name', 'original_taxon_name', 'iso', 'concept', 'ipa', 'tokens', 'alignment',
       'original_alignment', 'cogid', 'note']

wl = Wordlist(D)

# iterate over alignments and adjust them
for concept in wl.concepts:
    idxs = wl.get_list(concept=concept, flat=True)
    
    alms = [wl[idx,'alignment'].split(' ') for idx in idxs]
    max_rows = max([len(row) for row in alms])
    for i,row in enumerate(alms):
        if len(row) != max_rows:
            alms[i] = row + ['-' for x in range(max_rows)]
            alms[i] = alms[i][:max_rows]
            
    calms = misc.transpose(alms)

    out = []
    for i,row in enumerate(calms):
        if len(row) == row.count('-'):
            out += [i]
        elif '(' in row and '-' in row and len(set(row)) == 2:
            calms[i] = ['(' for j in range(len(alms))]
        elif ')' in row and '-' in row and len(set(row)) == 2:
            calms[i] = [')' for j in range(len(alms))]
    
    for i in out[::-1]:
        del calms[i]

    nalms = misc.transpose(calms)
    
    for i,idx in enumerate(idxs):
        wl[idx][wl.header['alignment']] = ' '.join(nalms[i])
        
alm = Alignments(wl)
alm.output('tsv', filename='payne-maipuran', ignore='all', prettify=False)

