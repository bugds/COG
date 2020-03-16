accessionNumbers = 'merged.txt'
gene2refseq = '/windows1/usr/Boog/gene2refseq/gene2refseq'
taxid2name = '/windows1/usr/Boog/taxid2name/s-names.dmp'
treeIn = 'mafft_fasta_treefile.treefile'
treeOut = '.'.join(treeIn.split('.')[:-1]) + '_reformat.nwk'

with open(accessionNumbers, 'r') as inputFile, \
    open(gene2refseq, 'r') as g2r, \
    open(taxid2name, 'r') as t2n, \
    open(treeIn, 'r') as nwk, \
    open(treeOut, 'w') as outFile:
    pDict = dict()
    line = inputFile.readline().replace('\n', '')
    while line:
        pDict[line] = ''
        line = inputFile.readline().replace('\n', '')
    
    species = dict()
    
    line = t2n.readline()
    while line:
        species[line.split('\t')[0]] = line.split('\t')[-1]
        line = t2n.readline()
    
    line = g2r.readline()
    while line:
        if line.split('\t')[5] in pDict:
            pDict[line.split('\t')[5]] = '{}_-_{}_{}'.format(
                line.split('\t')[5],
                species[line.split('\t')[0]].replace(' ', '_'),
                line.split('\t')[-1]
            )
        line = g2r.readline()
    
    tree = nwk.read()
    
    for key, value in pDict.items():
        tree = tree.replace(key, value)
    
    outFile.write(tree)
