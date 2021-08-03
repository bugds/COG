import sys
import os
import pickle
import subprocess
import datetime
import networkx
from io import StringIO
from Bio import SearchIO
from copy import deepcopy
from networkx.algorithms import clique
from pyvis.network import Network

# 'rootFolder' is a directory that contains:
# /preInput             (accession numbers of queried proteins, each in a 
#                       separate file named accordingly)
# /Input                (hits of initial Blast search)
# /Previous_Proteins    (dictionary of candidate proteins with metadata)
# /Blast_XML            (search results in xml format)
# /Results              (for output)
rootFolder = sys.path[0]
# ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
path2G2R = '/home/bioinfuser/data/corgi_files/corgi_oct/gene2refseq' 
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
path2T2N = '/home/bioinfuser/data/corgi_files/corgi_oct/names.dmp'
# Name of database with representative taxids
databaseName = 'clcn'
# Path to Blastp utility
path2blastp = '/home/bioinfuser/applications/ncbi-blast-2.10.1+/bin/blastp'

# E-value is generally 10^-4
evalueLimit = 0.1
# Query cover: length of a domain of interest divided by query length
qCoverLimit = 0.1
# Number of targets in initial Blast search (expected number of homologs)
initBlastTargets = '700'
# Number of CPU threads
numThreads = '48'
# Technical constant, do not change
blastChunkSize = 100
# A fraction of isoforms needed to be the closest between two genes,
# so the genes can be called homologous
orthologyThreshold = 1.0
# First step: initial Blast search, creating a dictionary of candidate-proteins
preInput = False
# Second step: merging results of Blast search (optional)
mergeInput = False
# Third step: perform Blast search, create dictionary of results
doMainAnalysis = False
# Forth step: analysis
finalAnalysis = True
# Remove all Blast results (to save space)
removeXml = True

class ProteinClass():
    '''Class for proteins
    '''
    def __init__(self, species, taxid, symbol, gene, refseq):
        '''Initialization
        :param species: Species in which proteins are synthetized
        :param taxid: Taxid of a species
        :param symbol: Gene symbol
        :param gene: Coding gene
        :param refseq: Reference sequence accession number
        '''
        self.species = species
        self.taxid = taxid
        self.symbol = symbol
        self.gene = gene
        self.refseq = refseq
        self.good = False # This parameter defines orthologs

def createInputForBlast(inp, filename):
    '''Creating file for Blastp input
    :param inp: Contents of the input
    :param filename: Name of currently analyzed file
    :return: Path to the temporary file
    '''
    with open('{}/Temp/{}{}'.format(rootFolder, filename, '.q'), 'w') as f:
        if isinstance(inp, str):
            f.write(inp)
        else:
            f.write('\n'.join(inp))
    return '{}/Temp/{}{}'.format(rootFolder, filename, '.q')

def bashBlast(
    query, 
    out, 
    outfmt='5', 
    num_threads=numThreads, 
    max_target_seqs='500'):
    '''Run Blastp search
    :param query: File with query accession number
    :param out: Output file path
    :param outfmt: Output format
    :param num_threads: Number of threads
    :max_target_seqs: Number of target sequences
    :return: False if failed, True if succeded
    '''
    blastProcess = subprocess.run(
        [path2blastp,
        '-db', databaseName,
        '-query', query,
        '-outfmt', '5',
        '-out', out,
        '-num_threads', num_threads,
        '-max_target_seqs', max_target_seqs],
        stderr = subprocess.PIPE
    )
    if 'Sequence ID not found' in blastProcess.stderr.decode():
        print(blastProcess.stderr.decode())
        exit(1)
        return False
    return True

def initialBlast(filename, query):
    '''Run initial Blast - results will constitute the list of
    candidate-homologs
    :param filename: Name of analyzed file
    :param query: Accession number of a query
    :return: Blast search results in xml-format
    '''
    query = createInputForBlast(query, filename)
    xmlPath = rootFolder + '/Blast_XML/' + os.path.splitext(filename)[0] + '.xml'
    bashBlast(
        query=query, 
        out=xmlPath,
        max_target_seqs=initBlastTargets
    )
    return SearchIO.parse(xmlPath, 'blast-xml')

def parseInitialBlast(blast):
    '''Filter results based on query cover and E-value limits
    :param blast: Blast search results in xml-format
    :return: Filtered list of proteins
    '''
    initBlastList = list()
    for record in blast:
        queryLen = int(record.seq_len)
        for hit in record:
            for hsp in hit:
                alnSpan = int(hsp.query_span)
                qCover = float("{0:.2f}".format(alnSpan/queryLen))
                if (qCover > qCoverLimit) and (hsp.evalue < evalueLimit):
                    # can't just take hit.accession - 
                    # does not have accession version
                    substrings = hit.id.split('|')
                    for i in range(len(substrings)):
                        if substrings[i] == 'ref':
                            initBlastList.append(substrings[i+1])
    return initBlastList

def checkPreviousPickle(filename, folder):
    '''Take previously pickled objects or returns "False"
    :param filename: Name of analyzed file
    :param folder: Folder in which pickled objects are contained
    :return: Object or "False" if it does not exist
    '''
    for prevName in os.listdir(rootFolder + folder):
        basePrevName = os.path.splitext(prevName)[0]
        baseFilename = os.path.splitext(filename)[0]
        if baseFilename == basePrevName:
            path = rootFolder + folder + '/' + prevName
            with open(path, 'rb') as f:
                return pickle.load(f)
    return False

def savePickle(shortName, toSave, folder):
    '''Save variables into a pickle file
    :param shortName: Part of the analyzed file name
    :param toSave: Object to save
    :param folder: Folder in which pickled objects are contained
    '''
    path = rootFolder + folder + '/' + shortName + '.pkl'
    with open(path, 'wb') as f:
        pickle.dump(toSave, f)

def getSequences(seqFilename, proteins):
    '''Get accession numbers from corresponding file
    :param seqFilename: Name of a file with accession numbers
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with proteins metadata
    '''
    path = rootFolder + '/Input/' + seqFilename
    seqFile = open(path, 'r')
    line = seqFile.readline().replace('\n', '')
    while line:
        if not line in proteins:
            proteins[line] = ProteinClass(None, None, None, None, None)
        line = seqFile.readline().replace('\n', '')
    return proteins

def parseG2RHeader(header):
    '''Get the columns of taxids, gene symbols, accession numbers, and
    gene ID from gene2refseq file
    :param header: Header of a gene2refseq file
    :return: Dictionary of column numbers
    '''
    return {
            't':header.index('#tax_id'), 
            'g':header.index('GeneID') , 
            'p':header.index('protein_accession.version') , 
            's':header.index('Symbol')
        }

def saveTempSet(toSave, tempSet, proteins, index):
    '''Save set of data for a single protein, parsed from
    gene2refseq file
    :param toSave: If False, function will not be performed
    :param tempSet: Data for a currently parsed gene
    :param proteins: Dictionary for storing information about proteins
    :param index: Dictionary of gene2refseq columns
    :return: Supplemented dictionary for storing information about proteins
    '''
    if toSave:
        for l in tempSet:
            if l.split('\t')[index['p']] != '-':
                proteins[l.split('\t')[index['p']]] = ProteinClass(
                    None,
                    l.split('\t')[index['t']], 
                    l.split('\t')[index['s']],
                    l.split('\t')[index['g']], 
                    l.split('\t')[index['p']]
                )
    return proteins

def getIsoforms(proteins):
    '''Get isoforms for all featured proteins
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with isoforms
    '''
    with open(path2G2R, 'r') as gene2Refseq:
        tempSet = [gene2Refseq.readline().replace('\n', '')]
        index = parseG2RHeader(tempSet[0].split('\t'))
        line = gene2Refseq.readline().replace('\n', '')
        toSave = False
        while line:
            if tempSet[-1].split('\t')[index['p']] in proteins:
                toSave = True
            if line.split('\t')[index['g']] == tempSet[-1].split('\t')[index['g']]:
                tempSet.append(line)
            else:
                saveTempSet(toSave, tempSet, proteins, index)
                toSave = False
                tempSet = [line]
            line = gene2Refseq.readline().replace('\n', '')
        saveTempSet(toSave, tempSet, proteins, index)
    return proteins

def getSpeciesName(proteins):
    '''Get names of species from taxids, using the names.dmp file
    :param proteins: Dictionary for storing information about proteins
    :return: Supplemented dictionary for storing information about proteins
    '''
    with open(path2T2N, 'r') as f:
        taxids = [p.taxid for p in proteins.values()]
        line = f.readline()
        while line:
            lineList = line.split('\t')
            if lineList[0] in taxids:
                if 'scientific name' in lineList[6]:
                    for p in proteins.values():
                        if p.taxid == lineList[0]:
                            p.species = lineList[2]
                            if 'gorilla' in p.species:
                                p.species = 'Gorilla gorilla gorilla'
            line = f.readline()
    return proteins

def writeInBlastDict(blast, blastDict):
    '''Create or append to a dictionary containing Blast results
    :param blast: Contents of the XML-file with Blast results
    :param blastDict: Dictionary containing part of Blast results (or empty)
    :return: Dictionary containing Blast results
    '''
    for record in blast:
        if record.id not in blastDict:
            blastDict[record.id] = {}
        for hit in record:
            species = hit.description.split('[')[1].split(']')[0]
            if not species in blastDict[record.id]:
                substrings = hit.id.split('|')
                for i in range(len(substrings)):
                    if substrings[i] == 'ref':
                        blastDict[record.id][species] = substrings[i+1]
    return blastDict
    
def blastSearch(query, speciesList, filename, blastDict):
    '''Run Blast, save results of a search to a file and return its contents
    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which Blast is performed
    :param filename: Name of original fasta file for saving results of Blast
    :return: Contents of xml-file with Blast results in form of a dictionary
    '''
    xmlPath = rootFolder \
        + '/Blast_XML/' \
        + os.path.splitext(filename)[0] \
        + '.xml'
    query = createInputForBlast(query, filename)
    blastNotVoid = bashBlast(
        query=query,
        out=xmlPath,
    )
    if blastNotVoid:
        blast = SearchIO.parse(xmlPath, 'blast-xml')
        writeInBlastDict(blast, blastDict)
    os.remove(query)
    if removeXml:
        os.remove(xmlPath)
    return blastDict

def clearProteins(proteins):
    '''Delete proteins for which metadata was not found
    :param proteins: Dictionary for storing information about proteins
    :return: Shortened proteins dictionary
    '''
    toDel = list()
    for r in proteins.keys():
        if proteins[r].species == None:
            toDel.append(r)
            print('No info on protein ' + r)
    for r in toDel:
        del proteins[r]
    return proteins

def createBlastDict(proteins, filename):
    '''Run Blast algorithm to create dictionary based on its results
    :param proteins: Dictionary for storing information about proteins
    :param filename: Name of currently analyzed file
    :return: Dictionary containing Blast results
    '''
    blastDict = dict()
    chunksForBlast = dict()
    chunkN = 0
    counter = 0
    for p in proteins.values():
        chunksForBlast[p.refseq] = p
        counter += 1
        if counter >= blastChunkSize:
            blastDict = blastSearch(
                sorted([seq.refseq for seq in chunksForBlast.values()]),
                sorted([seq.taxid for seq in proteins.values()]),
                '{}_basic_chunk{}'.format(
                    os.path.splitext(filename)[0],
                    str(chunkN) + '.nomatter'
                ),
                blastDict
            )
            print(
                str(datetime.datetime.now()) \
                + ': Blast search completed (chunk ' \
                + str(chunkN) \
                + ')' \
            )
            counter = 0
            chunkN += 1
            chunksForBlast = dict()
            savePickle('part_' + os.path.splitext(filename)[0], \
                {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
    blastDict = blastSearch(
        sorted([seq.refseq for seq in chunksForBlast.values()]),
        sorted([seq.taxid for seq in proteins.values()]),
        '{}_basic_chunk{}'.format(
            os.path.splitext(filename)[0],
            str(chunkN) + '.nomatter'
        ),
        blastDict
    )
    print(
        str(datetime.datetime.now()) \
        + ': Blast search completed (chunk ' \
        + str(chunkN) \
        + ')' \
    )
    savePickle(os.path.splitext(filename)[0], \
        {'proteins':proteins, 'blastDict':blastDict}, '/For_online')

    print('Checking Blast dictionary...') 
    blastDict = checkBlastDict(proteins, filename, blastDict, 0)
    print(str(datetime.datetime.now()) + ': Blast dictionary checked')
    savePickle(os.path.splitext(filename)[0], \
        {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
    return blastDict

def checkBlastDict(proteins, filename, blastDict, iteration, previous=[set(), set()]):
    '''Check if Blast found all species in each case
    :param filename: Name of currently analyzed file
    :param blastDict: Dictionary containing Blast results
    :param proteins: Dictionary for storing information about proteins
    :param iteration: Number of additional Blast required currently
    :return: "blastDict" supported with Blast results
    '''
    queriesForBlast = set()
    taxidsForBlast = set()
    for seq in proteins.keys():
        taxidsAll = set([p.taxid for p in proteins.values()])
        taxidsAlreadyIn = set([
            p.taxid for p in proteins.values() \
            if p.species in blastDict[seq]
        ])
        if (taxidsAll - taxidsAlreadyIn):
            queriesForBlast.add(seq)
            taxidsForBlast = taxidsForBlast | (taxidsAll - taxidsAlreadyIn)
    if (previous[0] == queriesForBlast) and (previous[1] == taxidsForBlast):
        for q in queriesForBlast:
            for t in taxidsForBlast:
                s = [p.species for p in proteins.values() if p.taxid == t][0]
                if not s in blastDict[q]:
                    blastDict[q][s] = 'NA'
        return blastDict
    else:
        chunksForBlast = dict()
        counter = 0
        chunkN = 0
        for refseq in queriesForBlast:
            chunksForBlast[refseq] = proteins[refseq]
            counter += 1
            if counter >= blastChunkSize:
                blastDict = blastSearch(
                    sorted(list([seq.refseq for seq in chunksForBlast.values()])),
                    sorted(list(taxidsForBlast)),
                    '{}_check_chunk{}_iter{}'.format(
                        os.path.splitext(filename)[0],
                        str(chunkN),
                        str(iteration) + '.nomatter'
                    ),
                    blastDict
                )
                print(
                    str(datetime.datetime.now()) \
                    + ': Blast search completed (iteration ' \
                    + str(iteration) \
                    + ', chunk ' \
                    + str(chunkN) \
                    + ')' \
                )
                counter = 0
                chunkN += 1
                chunksForBlast = dict()
                savePickle('part_' + os.path.splitext(filename)[0], \
                    {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
        blastDict = blastSearch(
            sorted(list([seq.refseq for seq in chunksForBlast.values()])),
            sorted(list(taxidsForBlast)),
            '{}_check_chunk{}_iter{}'.format(
                os.path.splitext(filename)[0],
                str(chunkN),
                str(iteration) + '.nomatter'
            ),
            blastDict
        )
        print(
            str(datetime.datetime.now()) \
            + ': Blast search completed (iteration ' \
            + str(iteration) \
            + ', chunk ' \
            + str(chunkN) \
            + ')' \
        )
        savePickle(os.path.splitext(filename)[0], \
            {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
        return checkBlastDict(proteins, filename, blastDict, iteration + 1,
        [queriesForBlast, taxidsForBlast])

def createDictsForAnalysis(proteins, blastDict):
    '''Create a deep copy of blastDict for modifications, as well as
    dictionary of genes, analogous to blastDict (which is for proteins),
    applying the "orthologyThreshold" variable
    :param blastDict: Dictionary containing Blast results
    :param proteins: Dictionary for storing information about proteins
    :return: Deep copy of blastDict and blastDict for genes
    '''
    transDict = deepcopy(blastDict)
    for q in transDict.keys():
        for s in transDict[q].keys():
            if transDict[q][s] in proteins:
                transDict[q][s] = proteins[transDict[q][s]].gene
            else:
                transDict[q][s] = 'NA'
    geneDict = dict()
    for g in set([p.gene for p in proteins.values()]):
        geneDict[g] = dict()
        isoforms = [p.refseq for p in proteins.values() if p.gene == g]
        for s in set([p.species for p in proteins.values()]):
            targetGenes = dict()
            for i in isoforms:
                if s in transDict[i]:
                    if not transDict[i][s] in geneDict[g]:
                        targetGenes[transDict[i][s]] = 1
                    else:
                        targetGenes[transDict[i][s]] += 1
            if len(targetGenes) > 0:
                if max(targetGenes.values())/sum(targetGenes.values()) >= orthologyThreshold:
                    maxKeys = [k for k, v in targetGenes.items() if v == max(targetGenes.values())]
                    if len(maxKeys) != 1:
                        print('Multiple equally good orthologous genes for ' + g + ': ' + ', '.join(maxKeys) + '. Chosen ' + maxKeys[0])
                    geneDict[g][s] = maxKeys[0]
    return transDict, geneDict

def findLargestMaxCliques(graph, mainGene):
    '''Find largest maximal cliques
    :param graph: Graph of Blast results
    :param mainGene: Gene of query
    :return: Largest maximal cliques list
    '''
    maxLen = clique.node_clique_number(graph, mainGene)
    maxCliques = list()
    for c in clique.find_cliques(graph):
        if len(c) == maxLen:
            if mainGene in c:
                maxCliques.append(c)
    return maxCliques

def createGraph(mainGene, mainSpecies, proteins, geneDict):
    '''Create graph
    :param mainGene: Gene of query
    :param mainSpecies: Species of query
    :param proteins: Dictionary for storing information about proteins
    :param geneDict: Dictionary of genes according to Blast results
    :return: Graph representating Blast results
    '''
#    graph = networkx.Graph()
#    graph.add_node(mainGene)
#    for q in geneDict:
#        qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
#        if (mainSpecies in geneDict[q]) and (qSpecies in geneDict[mainGene]):
#            if (geneDict[q][mainSpecies] == mainGene) and \
#            (geneDict[mainGene][qSpecies] == q):
#                graph.add_node(q)
#    for q in graph.nodes():
#        for t in graph.nodes():
#            qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
#            tSpecies = [p.species for p in proteins.values() if p.gene == t][0]
#            if (tSpecies in geneDict[q]) and (qSpecies in geneDict[t]):
#                if (q != t) and (geneDict[q][tSpecies] == t) and (geneDict[t][qSpecies] == q):
#                    graph.add_edge(q, t)
#    maxCliques = findLargestMaxCliques(graph, mainGene)
#    return graph, maxCliques
    graph = networkx.Graph()
    graph.add_node(mainGene)
    for q in geneDict:
        graph.add_node(q)
    for q in graph.nodes():
        for t in graph.nodes():
            qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
            tSpecies = [p.species for p in proteins.values() if p.gene == t][0]
            if (tSpecies in geneDict[q]) and (qSpecies in geneDict[t]):
                if (q != t) and (geneDict[q][tSpecies] == t) and (geneDict[t][qSpecies] == q):
                    graph.add_edge(q, t)
    maxCliques = findLargestMaxCliques(graph, mainGene)
    return graph, maxCliques

def drawGraph(
    graph, 
    maxCliques,
    proteins, 
    filename, 
    mainGene,
    mainSpecies,
    springLength = 200,
    common_color = 'rgb(30 ,144 ,255)',
    main_color = 'red',
    max_color = 'rgb(50,205,50)'):
    '''Draw graph
    :param graph: Graph representing Blast results
    :param proteins: Dictionary for storing information about proteins
    :param filename: Name of analyzed file
    :param mainGene: Gene of query
    '''
    net = Network(height = '95%', width = '100%')
    net.from_nx(graph)
    net.barnes_hut()
    net.repulsion(spring_length = springLength)
#    maxLen = clique.node_clique_number(graph, mainGene)
#    cliquesNodes = {n:clique.node_clique_number(graph, n) for n in graph.nodes()}
#    maxNodWes = [n for n in cliquesNodes.keys() if cliquesNodes[n] == maxLen]
    moreMaxCliques = [g for c in maxCliques for g in c]
    for G in [p.gene for p in proteins.values() if p.species == mainSpecies]:
        moreMaxCliques.append([g for c in findLargestMaxCliques(graph, G) for g in c])
    maxNodes = list(set().union(*moreMaxCliques))
    for node in net.nodes:
        node['group'] = 'common'
        node['color'] = common_color
#        neighbors = [n for n in graph.neighbors(node['label'])]
#        percentage = int(100*len(neighbors)/(len(graph) - 1))
#        node['title'] = 'Connected to ' + str(percentage) + '% nodes.'
        node['title'] = node['label']
        node['size'] = 3
        node['font'] = dict()
        node['font']['size'] = 8
    for node in net.nodes:
        if node['label'] \
          in [p.gene for p in proteins.values() if p.species == mainSpecies]:
            node['group'] = 'main'
            node['color'] = main_color
        elif node['label'] in maxNodes:
            node['group'] = 'max'
            node['color'] = max_color
        node['label'] = '_'.join([\
            [p.species for p in proteins.values() if p.gene == node['label']][0],
            [p.symbol for p in proteins.values() if p.gene == node['label']][0]
        ]).replace(' ', '_')
    for edge in net.edges:
        edge['width'] = 0
        edge['color'] = dict()
        edge['color'] = common_color
        if edge['from'] in maxNodes:
            edge['color'] = max_color
        if edge['to'] in maxNodes:
            edge['color'] = max_color
        if edge['from'] in \
          [p.gene for p in proteins.values() if p.species == mainSpecies]:
            edge['color'] = main_color
        if edge['to'] in \
          [p.gene for p in proteins.values() if p.species == mainSpecies]:
            edge['color'] = main_color
    net.save_graph(rootFolder + '/Results/' + os.path.splitext(filename)[0] + '_pyvis.html')

def goodGeneMakesGoodProtein(proteins, goodGenes):
    '''If any isoform of a gene hypothesized to be orthologous, all
    other isoforms of this gene should be considered orthologous
    :param proteins: Dictionary for storing information about proteins
    :param goodGenes: Genes, isoforms of which are hypothesized to be orthologous
    :return: Changed "proteins" dictionary
    '''
    for g in goodGenes:
        isoforms = [p.refseq for p in proteins.values() if p.gene == g]
        for i in isoforms:
            proteins[i].good = True
    return proteins

def clearProteins2(proteins):
    '''Create a deep copy of proteins and clear it from proteins that
    are considered non-orthologous
    :param proteins: Dictionary for storing information about proteins
    :return: Changed copy of "proteins" dictionary
    '''
    refDict = dict()
    for p in proteins.values():
        if p.good:
            if p.species in refDict:
                if refDict[p.species] != p.gene:
                    print('Multiple genes of ' + p.species + ' are good')
            refDict[p.species] = p.gene                
    toDel = set()
    for p in proteins.values():
        if (p.species in refDict.keys()) and (p.gene != refDict[p.species]):
            toDel.add(p.refseq)
    tempProteins = deepcopy(proteins)
    for refseq in toDel:
        tempProteins.pop(refseq)
    gSpecies = set([tempProteins[k].species for k in tempProteins if tempProteins[k].good])
    bSpecies = set([tempProteins[k].species for k in tempProteins if not tempProteins[k].good])
    if bSpecies.intersection(gSpecies) != set():
        raise Exception('Same species in referencial and non-referencial groups!')
    return tempProteins

def writeHtml(
    proteins,
    gSpecies,
    gRefseqs,
    qSpecies,
    qGenes,
    qReverse,
    qForward):
    '''Write Blast analysis for single gene in form of an HTML-file
    :param proteins: Dictionary for storing information about proteins
    :param gSpecies: Set of good species
    :param gRefseqs: Set of good accession numbers
    :param qSpecies: Species linked to analyzed genes
    :param qGenes: Set of analyzed genes
    :param qReverse: Dictionary of reciprocal Blast: isoform>species
    :param qForward: Set of species found by forward Blast
    :return: HTML-string of Blast analysis for single species
    '''
    htmlPart = StringIO()
    htmlString = list()
    htmlString.append('<details>\n\t<summary>{}</summary>\n')
    htmlString.append('\t<details>\n\t\t<summary>&emsp;Gene id: {}</summary>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{}/{} referencial -> this gene. Fails:</summary>\n')
    htmlString.append('\t\t\t\t&emsp;&emsp;&emsp;&emsp;{} [{}]<br>\n')
    htmlString.append('\t\t</details>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{}/{} isoforms of this gene -> only referencial isoforms. Fails:</summary>\n')
    htmlString.append('\t\t\t<details>\n\t\t\t\t<summary>&emsp;&emsp;&emsp;{} -> isoforms of {}/{} referencial genes. Fails:</summary>\n')
    htmlString.append('\t\t\t\t\t&emsp;&emsp;&emsp;&emsp;&emsp;{}<br>\n')
    htmlString.append('\t\t\t</details>\n')
    htmlString.append('\t\t</details>\n\t</details>\n')
    htmlString.append('</details>')
    htmlPart.write(htmlString[0].format(qSpecies))
    for qGene in sorted(list(qGenes)):
        temporaryProteins = [p for p in proteins if proteins[p].gene == qGene]
        htmlPart.write(htmlString[1].format(
        ', Gene symbol: '.join([qGene, proteins[temporaryProteins[0]].symbol]),
        len(qForward[qGene]),
        len(gRefseqs)
        ))
        for fail in sorted(list((gRefseqs - qForward[qGene]))):
            htmlPart.write(htmlString[2].format(
                fail,
                proteins[fail].species
            ))
        htmlPart.write(htmlString[3].format(
            len([qR for qR in qReverse[qGene].values() if len(qR) == len(gSpecies)]),
            len(qReverse[qGene])
        ))
        for isoform, success in qReverse[qGene].items():
            htmlPart.write(htmlString[4].format(
                isoform,
                len(success),
                len(gSpecies)
            ))
            for fail in sorted(list((gSpecies - success))):
                htmlPart.write(htmlString[5].format(
                    fail
                ))
            htmlPart.write(htmlString[6])
        htmlPart.write(htmlString[7])
    htmlPart.write(htmlString[8])
    return htmlPart.getvalue()

def analyzeBlastDict(blastDict, proteins):
    '''Analysis of a Blast dictionary
    :param blastDict: Dictionary containing Blast results
    :param proteins: Dictionary for storing information about proteins
    :return: HTML-string containing analysis results
    '''
    htmlFull = ''
    gProteins = [p for p in proteins.values() if p.good]
    gSpecies = set([p.species for p in gProteins])
    for qSpecies in sorted(list(set([p.species for p in proteins.values()]))):
        qGenes = set()
        qReverse = dict()
        qForward = dict()
        for qGene in sorted([p.gene for p in proteins.values() if p.species == qSpecies]):
            qGenes.add(qGene)
            qReverse[qGene] = dict() # for qRefseq use qReverse.keys()
            qForward[qGene] = set()
            # Reciprocal Blast
            for qRefseq in sorted([p.refseq for p in proteins.values() if p.gene == qGene]):
                #if not proteins[qRefseq].good:
                qReverse[qGene][qRefseq] = set()
                for s in sorted(gSpecies):
                    if blastDict[qRefseq][s] in [p.refseq for p in gProteins]:
                        qReverse[qGene][qRefseq].add(s)
            # Forward Blast
            for gRefseq in sorted([p.refseq for p in gProteins]):
                if blastDict[gRefseq][qSpecies] in proteins:
                    if proteins[blastDict[gRefseq][qSpecies]].gene == qGene:
                        qForward[qGene].add(gRefseq)
        htmlFull += writeHtml(
                proteins,
                gSpecies,
                set([p.refseq for p in gProteins]),
                qSpecies,
                qGenes,
                qReverse,
                qForward
        )
    return htmlFull

def reportHtml(filename, maxClique, proteins, html):
    '''Write results in form of HTML-file
    :param filename: Name of analyzed file
    :param maxClique: Currently analyzed largest maximal clique
    :param proteins: Dictionary for storing information about proteins
    :param html: A string containing the results of an algorithm
    '''
    with open(
        rootFolder \
        + '/Results/' \
        + os.path.splitext(filename)[0] \
        + '.html', 
        'w'
    ) as out:
        out.write('Original cluster:<br>')
        for gene in sorted(list(maxClique)):
            out.write(
                [p.species for p in proteins.values() \
                    if p.gene == gene][0] + ': ' + \
                [p.symbol for p in proteins.values() \
                    if p.gene == gene][0] + ' (' + \
                gene + ')<br>')
        out.write('<br>Results:<br>')
        out.write(html)

def runPreInput():
    '''Run the first step -
    query Blast, get orthologs-candidates
    '''
    for filename in os.listdir(rootFolder + '/preInput'):
        print(filename)
        with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
            mainRefseq = oneStrFile.read().replace('\n', '')
            blast = initialBlast(filename, mainRefseq)
            initBlastList = parseInitialBlast(blast)
            with open(rootFolder + '/Input/' + filename, 'w') as blastResults:
                blastResults.write('\n'.join(list(dict.fromkeys(initBlastList))))

def runMergeInput():
    '''Run the second step (optional) -
    merge all results into one protein database
    '''
    mergedSet = set()
    for filename in os.listdir(rootFolder + '/Input'):
        with open(rootFolder + '/Input/' + filename, 'r') as singleFile:
            singleContent = singleFile.read()
        mergedSet = mergedSet | set(singleContent.split('\n'))
    mergedSet.discard('')
    with open(rootFolder + '/Input/merged.txt', 'w') as mergedFile:
        mergedFile.write('\n'.join(mergedSet))

def runMainAnalysis():
    '''Run the third step -
    enrichment of protein database and Blast search
    '''
    for filename in os.listdir(rootFolder + '/Input'):
        proteins = getSequences(filename, dict())
        proteins = getIsoforms(proteins)
        proteins = getSpeciesName(proteins)
        proteins = clearProteins(proteins)
        savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')
        print(str(datetime.datetime.now()) + ': "proteins" ready')
        blastDict = createBlastDict(proteins, filename)

def runFinalAnalysis():
    '''Run the forth step -
    analysis of Blast results
    '''
    for filename in os.listdir(rootFolder + '/preInput'):
        print(filename)
        # proteins need to be refreshed each time we do an analysis
        # else good values are not dropped
        pkl = checkPreviousPickle(filename, '/For_online')
        blastDict = pkl['blastDict']
        pkl = checkPreviousPickle(filename, '/Previous_Proteins')
        proteins = pkl
        transDict, geneDict = createDictsForAnalysis(proteins, blastDict)
        with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
            mainRefseq = oneStrFile.read().replace('\n', '')
        for k in proteins.keys():
            if k.split('.')[0] == mainRefseq:
                mainRefseq = k
        mainSpecies = proteins[mainRefseq].species
        mainGene = proteins[mainRefseq].gene 
        graph, maxCliques = createGraph(mainGene, mainSpecies, proteins, geneDict)
        drawGraph(graph, maxCliques, proteins, filename, mainGene, mainSpecies)
        maxCliques.sort()
        cliqueCounter = 0
        for maxClique in maxCliques:
            cliqueCounter += 1
            for p in proteins:
                proteins[p].good = False
            proteins = goodGeneMakesGoodProtein(proteins, maxClique)
            tempProteins = clearProteins2(proteins)
            html = analyzeBlastDict(blastDict, tempProteins)
            reportHtml(os.path.splitext(filename)[0] + '_clique' + str(cliqueCounter), maxClique, proteins, html)

def main():
    '''Main function
    '''
    print(str(datetime.datetime.now()) + ': start')
    if preInput:
        runPreInput()
    if mergeInput:
        runMergeInput()
    if doMainAnalysis:
        runMainAnalysis()
    if finalAnalysis:
        runFinalAnalysis()

main()
