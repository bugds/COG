#! /usr/bin/env python3

import sys
import os
import json
import pickle
import subprocess
import datetime
import networkx
import argparse
import markov_clustering
from io import StringIO
from Bio import SearchIO
from copy import deepcopy
from random import choice
from networkx.algorithms import clique
from pyvis.network import Network

parser = argparse.ArgumentParser()
parser.add_argument(
    'input', type=str,
    help = "Input directory (separate files with RefSeq identifier in each)")
parser.add_argument(
    'output', type=str,
    help = "Output directory")
parser.add_argument(
    "-e", "--eval", type = float,
    help = "E-value limit", default = 0.0001)
parser.add_argument(
    "-q", "--qcov", type = float,
    help = "Query cover limit", default = 0.1)
parser.add_argument(
    "-b", "--init", type = str,
    help = "Number of initial BLAST targets", default = '1500')
parser.add_argument(
    "-t", "--threadNum", type = str,
    help = "Number of threads", default = '40')
parser.add_argument(
    "-o", "--orthology", type = float,
    help = "Orthology threshold", default = 1)
parser.add_argument(
    "-s", "--stage", type = str,
    help = "Run a particular stage", default = 'all',
    choices = ["1", "2", "3", "all"])
parser.add_argument(
    "-m", "--merge",
    help = "Merge BLAST hits from multiple queries",
    default = False, action = 'store_true')
parser.add_argument(
    "-r", "--removeXML",
    help = "Don't remove XML files",
    default = True, action = 'store_false')
parser.add_argument(
    "-a", "--algorithm", type = str,
    help = "Graph building algorithm", default = 'strict',
    choices = ["best", "strict"])
parser.add_argument(
    "-g", "--gravity", type = float,
    help = "Gravitational constant", default = -30000)
parser.add_argument(
    "-c", "--config", type = str,
    help = "Configuration file", default = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "cogconf.txt"))
args = parser.parse_args()

with open(args.config, 'r') as inp:
    confLines = inp.readlines()
    for line in confLines:
        if line.startswith('path2G2R'):
            path2G2R = line.split(':')[1].strip()
        elif line.startswith('path2T2N'):
            path2T2N = line.split(':')[1].strip()
        elif line.startswith('databaseName'):
            databaseName = line.split(':')[1].strip()
        elif line.startswith('path2blastp'):
            path2blastp = line.split(':')[1].strip()
        elif line.startswith('blastdbcmd'):
            blastdbcmd = line.split(':')[1].strip()

# Input file
inputFile = args.input
# Output directory
rootFolder = args.output
# E-value is generally 10^-4
evalueLimit = float(args.eval)
# Query cover: length of a domain of interest divided by query length
qCoverLimit = float(args.qcov)
# Number of targets in initial Blast search (expected number of homologs)
initBlastTargets = str(args.init)
# Number of CPU threads
numThreads = str(args.threadNum)
# Technical constant, do not change
blastChunkSize = 100
# A fraction of isoforms needed to be the closest between two genes,
# so the genes can be called homologous
orthologyThreshold = float(args.orthology)
# Visualization constant
gravityConstant = float(args.gravity)

if str(args.stage) == "1":
    # First step: initial Blast search, creating a dictionary of candidate-proteins
    preInput = True
    doMainAnalysis = False
    finalAnalysis = False
elif str(args.stage) == "2":
    # Third step: perform Blast search, create dictionary of results
    preInput = False
    doMainAnalysis = True
    finalAnalysis = False
elif str(args.stage) == "3":
    # Forth step: analysis
    preInput = False
    doMainAnalysis = False
    finalAnalysis = True
else:
    preInput = True
    doMainAnalysis = True
    finalAnalysis = True

# Merging results of Blast search (optional)
mergeInput = args.merge
# Remove all Blast results (to save space)
removeXml = args.removeXML

# Defines if the best isoform or the majority of isoforms will indicate the orthology
algorithm = args.algorithm

os.makedirs(rootFolder, exist_ok = True)
inputDir = os.path.join(rootFolder, 'preInput')

if (len(os.listdir(rootFolder)) != 0):
    if (args.stage == "all") or (args.stage == "1"):
        raise Exception('Output directory not empty!')
elif (args.stage == "all") or (args.stage == "1"):
    os.makedirs(inputDir)
    os.makedirs(os.path.join(rootFolder, 'Blast_XML'))
    os.makedirs(os.path.join(rootFolder, 'For_online'))
    os.makedirs(os.path.join(rootFolder, 'Input'))
    os.makedirs(os.path.join(rootFolder, 'Previous_Proteins'))
    os.makedirs(os.path.join(rootFolder, 'Results'))
    os.makedirs(os.path.join(rootFolder, 'Temp'))

with open(inputFile, 'r') as inp:
    line = inp.readline().strip()
    while line:
        currPath = os.path.join(inputDir, line)
        if not(os.path.isfile(currPath)):
            with open(currPath, 'w') as out:
                out.write(line)
        else:
            with open(currPath, 'r') as prevInp:
                prevInputLines = prevInp.readlines()
                if line != prevInputLines[0].strip():
                    raise Exception('Wrong input file!')
                if len(prevInputLines) != 1:
                    raise Exception('Wrong input file!')
        line = inp.readline().strip()

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
        '-outfmt', outfmt,
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
    xmlPath = os.path.join(rootFolder, 'Blast_XML', os.path.splitext(filename)[0] + '.xml')
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
    for prevName in os.listdir(os.path.join(rootFolder, folder)):
        basePrevName = os.path.splitext(prevName)[0]
        baseFilename = os.path.splitext(filename)[0]
        if baseFilename == basePrevName:
            path = os.path.join(rootFolder, folder, prevName)
            with open(path, 'rb') as f:
                return pickle.load(f)
    return False

def savePickle(shortName, toSave, folder):
    '''Save variables into a pickle file
    :param shortName: Part of the analyzed file name
    :param toSave: Object to save
    :param folder: Folder in which pickled objects are contained
    '''
    path = os.path.join(rootFolder, folder, shortName + '.pkl')
    with open(path, 'wb') as f:
        pickle.dump(toSave, f)

def getSequences(seqFilename, proteins):
    '''Get accession numbers from corresponding file
    :param seqFilename: Name of a file with accession numbers
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with proteins metadata
    '''
    path = os.path.join(rootFolder, 'Input', seqFilename)
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
    xmlPath = os.path.join(rootFolder, 'Blast_XML', os.path.splitext(filename)[0] + '.xml')
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
                {'proteins':proteins, 'blastDict':blastDict}, 'For_online')
    if len(chunksForBlast) > 0:
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
        {'proteins':proteins, 'blastDict':blastDict}, 'For_online')

    print('Checking Blast dictionary...') 
    blastDict = checkBlastDict(proteins, filename, blastDict, 0)
    print(str(datetime.datetime.now()) + ': Blast dictionary checked')
    savePickle(os.path.splitext(filename)[0], \
        {'proteins':proteins, 'blastDict':blastDict}, 'For_online')
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
                    {'proteins':proteins, 'blastDict':blastDict}, 'For_online')
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
            {'proteins':proteins, 'blastDict':blastDict}, 'For_online')
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
    greatIso = dict()
    for q in blastDict.keys():
        for s in blastDict[q].keys():
            if not (s in greatIso):
                greatIso[s] = dict()
            h = blastDict[q][s]
            if h in proteins:
                g = proteins[h].gene
                if not (g in greatIso[s]):
                    greatIso[s][g] = dict()
                if not (h in greatIso[s][g]):
                    greatIso[s][g][h] = 0
                greatIso[s][g][h] += 1
    for s in greatIso:
        for g in greatIso[s]:
            currMax = 0
            for h in greatIso[s][g]:
                if greatIso[s][g][h] >= currMax:
                    currMax = greatIso[s][g][h]
                    currIso = h
            greatIso[s][g] = currIso
    geneDict = dict()
    for g in set([p.gene for p in proteins.values()]):
        geneDict[g] = dict()
        isoforms = [p.refseq for p in proteins.values() if p.gene == g]
        if algorithm == 'best':
            s = [p.species for p in proteins.values() if p.gene == g][0]
            if g in greatIso[s]:
                isoforms = [greatIso[s][g]]
            else:
                isoforms = []
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
                    maxKeys = [k for k in maxKeys if k != 'NA']
                    if len(maxKeys) > 1:
                        print('Multiple equally good orthologous genes for ' + g + ': ' + ', '.join(maxKeys) + '. Chosen ' + maxKeys[0])
                    elif len(maxKeys) == 0:
                        maxKeys = ['NA']
                    geneDict[g][s] = maxKeys[0]
    return transDict, geneDict, greatIso

def createFastasForTrees(proteins, greatIso, filename):
    with open(os.path.join(rootFolder, 'Temp', 'bdc.txt'), 'w') as f:
        for s in greatIso:
            for g in greatIso[s]:
                f.write(greatIso[s][g] + '\n')
    blastProcess = subprocess.run(
        [blastdbcmd,
        '-db', databaseName,
        '-entry_batch', os.path.join(rootFolder, 'Temp', 'bdc.txt'),
        '-out', os.path.join(rootFolder, 'Temp', 'bdc_out.txt')],
        stderr = subprocess.PIPE
    )
    with open(os.path.join(rootFolder, 'Temp', 'bdc_out.txt'), 'r') as f:
        fastaLines = f.readlines()
    with open(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '.fasta'), 'w') as o:
        for l in fastaLines:
            if l.startswith('>'):
                r = l.split()[0][1:]
                if len(proteins[r].symbol) > 0:
                    l = '>' + proteins[r].species + '_' + proteins[r].symbol + '\n'
                else:
                    l = '>' + proteins[r].species + '_' + 'GeneID' + str(proteins[r].gene) + '\n'
                l = l.replace(' ', '_')
            o.write(l)

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
    graph.remove_nodes_from(list(networkx.isolates(graph)))
    maxCliques = findLargestMaxCliques(graph, mainGene)
    return graph, maxCliques

def markovClustering(graph):
    '''Perform Markov clustering on a graph
    :param graph: Graph representing Blast results
    :return: Clusters
    '''
    matrix = networkx.to_scipy_sparse_array(graph)
    return markov_clustering.get_clusters(markov_clustering.run_mcl(matrix))

def drawGraph(
    graph, 
    maxCliques,
    proteins, 
    filename,
    mainSpecies,
    springLength = 100,
    commonColor = 'rgb(23,70,128)',
    mainColor = 'rgb(237,41,57)',
    maxColor = 'rgb(255,255,0)',
    voidColor = 'rgb(120,120,120)',
    palette = [
        'rgb(189,176,0)',
        'rgb(88,189,0)',
        'rgb(0,176,189)',
        'rgb(0,97,189)',
        'rgb(101,0,189)',
        'rgb(255, 236, 25)',
        'rgb(133, 255, 25)',
        'rgb(25, 240, 255)',
        'rgb(25, 144, 255)',
        'rgb(148, 25, 255)',
        'rgb(255, 25, 186)',
        'rgb(189,0,132)']):
    '''Draw graph
    :param graph: Graph representing Blast results
    :param proteins: Dictionary for storing information about proteins
    :param filename: Name of analyzed file
    :param mainSpecies: Species of query
    '''
    connectionsDict = dict()
    net = Network(height = '65%', width = '100%')
    clusters = markovClustering(graph)
    clusters.sort(key = len, reverse = True)
    nodesDict = {i: node for i, node in enumerate(graph.nodes())}
    net.from_nx(graph)
    net.barnes_hut(spring_length = springLength, gravity = gravityConstant)
    # net.show_buttons(filter_=['physics'])
    #net.repulsion(spring_length = springLength)
    moreMaxCliques = [g for c in maxCliques for g in c]
    for G in [p.gene for p in proteins.values() if p.species == mainSpecies]:
        if (G in graph):
            moreMaxCliques.append([g for c in findLargestMaxCliques(graph, G) for g in c])
    mainGenes = [p.gene for p in proteins.values() if p.species == mainSpecies]
    maxNodes = list(set().union(*moreMaxCliques))
    mainGenes = [p.gene for p in proteins.values() if p.species == mainSpecies]
    paletteDict = dict()
    paletteNum = 0
    displayedClustersNum = min(len(palette), len(clusters))
    for i in range(0, displayedClustersNum):
        for j in clusters[i]:
            if paletteNum < len(palette):
                paletteDict[nodesDict[j]] = palette[i]
    for node in net.nodes:
        node['title'] = node['label']
        node['physics'] = True
        node['hidden'] = False
        connectionsDict[node['title']] = 0
        for i in range(0, len(clusters)):
            for j in clusters[i]:
                if nodesDict[j] == node['title']:
                    node['markov'] = i
        if node['title'] in mainGenes:
            node['group'] = 'main'
            node['color'] = node['color0'] = node['color1'] = mainColor
        else:
            if node['title'] in paletteDict:
                node['color1'] = paletteDict[node['title']]
            else:
                node['color1'] = voidColor
            if node['title'] in maxNodes:
                node['group'] = 'max'
                node['color'] = node['color0'] = maxColor
            else:
                node['group'] = 'common'
                node['color'] = node['color0'] = commonColor
        node['label'] = '_'.join([\
            [p.species for p in proteins.values() if p.gene == node['label']][0],
            [p.symbol if len(p.symbol) > 0 else ('GeneID' + str(p.gene)) for p in proteins.values() if p.gene == node['label']][0]
        ]).replace(' ', '_')
    allNodes = set([n['title'] for n in net.nodes])
    for edge in net.edges:
        edge['width'] = 0
        if edge['to'] in maxNodes:
            edge['from'], edge['to'] = edge['to'], edge['from']
        if edge['to'] in mainGenes:
            edge['from'], edge['to'] = edge['to'], edge['from']
        for direction in ['from', 'to']:
            connectionsDict[edge[direction]] += 1
            allNodes.discard(edge[direction])
    for node in net.nodes:
        node['size'] = 2 + connectionsDict[node['title']] // 10
        node['borderWidth'] = 2 + connectionsDict[node['title']] // 40
        node['borderWidthSelected'] = 2 + connectionsDict[node['title']] // 10
        node['font'] = dict()
        node['font']['size'] = 3 + node['size']
        if node['title'] in allNodes:
            if node['color'] != mainColor:
                node['color'] = voidColor
    net.save_graph(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '_pyvis.html'))

def changeVisJS(filename):
    with open(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '_pyvis.html'), 'r') as f:
        htmlLines = f.readlines()
    with open(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '.fasta'), 'r') as f:
        fastaLines = f.readlines()
    fastaDict = dict()
    for l in fastaLines:
        if l.startswith('>'):
            currGene = l[1:].strip()
            fastaDict[currGene] = ''
        else:
            fastaDict[currGene] += l.strip()
    fastaJson = json.dumps(fastaDict)
    fastaJson = 'fasta = ' + fastaJson
    commentString = [htmlL for htmlL in htmlLines if 'from the python' in htmlL][0]
    nodesString = [htmlL for htmlL in htmlLines if 'nodes = new vis.DataSet' in htmlL][0]
    edgesString = [htmlL for htmlL in htmlLines if 'edges = new vis.DataSet' in htmlL][0]
    with open(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '_pyvis.html'), 'w') as f:
        f.write('<meta content="text/html;charset=utf-8" http-equiv="Content-Type">\n')
        f.write('<meta content="utf-8" http-equiv="encoding">\n')
        for l in htmlLines:
            # if ('<link rel="stylesheet" href=' in l) and ('vis' in l):
            #     l = l.split('href=')
            #     l[1] = '"https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis.css" type="text/css" />'
            #     l = 'href='.join(l)
            # if ('<script type="text/javascript" src=' in l) and ('vis' in l):
            #     l = l.split('scr=')
            #     l[1] = '"https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis-network.min.js"> </script>'
            #     l = 'scr='.join(l)
            if '#mynetwork {' in l:
                f.write('''
        /* The Modal (background) */
        .modal {
          display: none; /* Hidden by default */
          position: fixed; /* Stay in place */
          z-index: 1; /* Sit on top */
          padding-top: 100px; /* Location of the box */
          left: 0;
          top: 0;
          width: 100%; /* Full width */
          height: 100%; /* Full height */
          overflow: auto; /* Enable scroll if needed */
          background-color: rgb(0,0,0); /* Fallback color */
          background-color: rgba(0,0,0,0.4); /* Black w/ opacity */
        }
        
        /* Modal Content */
        .modal-content {
          background-color: #fefefe;
          margin: auto;
          padding: 20px;
          border: 1px solid #888;
          width: 80%;
          height: 70%;
          overflow-y: scroll;
        }

        .fasta {
          width: 100%;
          height: 90%;
          box-sizing: border-box;
          margin: 0px;
        }
        
        /* The Close Button */
        .close {
          color: #aaaaaa;
          right: 15%;
          font-size: 28px;
          font-weight: bold;
          position: fixed;
        }
        
        .close:hover,
        .close:focus {
          color: #000;
          text-decoration: none;
          cursor: pointer;
        }

''')
            if 'return network;' in l:
                f.write('''
        function clickEvent() {
            var options = {
                fields: ['id', 'font']
            };
            var currentlySelected = nodes.get(network.getSelectedNodes(), options);
            if (JSON.stringify(currentlySelected) != JSON.stringify(previouslySelected)) {
                for (var i = 0; i < previouslySelected.length; i++) {
                    previouslySelected[i].font.background = "transparent";
                }
                nodes.update(previouslySelected)
                for (var i = 0; i < currentlySelected.length; i++) {
                    currentlySelected[i].font.background = "rgb(220,220,220)";
                }
                nodes.update(currentlySelected)
                previouslySelected = currentlySelected;
            }
            addSelected();
        };
        
        network.on('dragStart', function() {
            clickEvent();
        });
        
        network.on('click', function() {
            clickEvent();
        });
''')
            if '"hideEdgesOnDrag": false,' in l:
                f.write('\t\t"multiselect": true,\n')
            if '"physics": {' in l:
                f.write('''
    "layout": {
        "improvedLayout": false
    },
''')
            if 'drawGraph();' in l:
                f.write('''
    // BRC code

    var labelIdDict = {};
    var clustalId = "";

    var ids = nodes.getIds();
    for (var i = 0; i < ids.length; i++) {
        labelIdDict[nodes.get(ids[i]).label] = ids[i];
    };
    
    function describe() {
        options = {
            fields: ['id', 'label', 'markov'],
            filter: function(item){
                return item.hidden == false;
            }
        };
        var editLabels = document.getElementById("editNodes").value.split("\\n");
        var editIds = [];
        for (var i = 0; i < editLabels.length; i++) {
            editIds.push(labelIdDict[editLabels[i]]);
        }
        editIds = editIds.filter(x => x !== undefined);
        if (editIds.length > 0) {
            editNodes = nodes.get(editIds, options)
            var description = '';
            for (var i = 0; i < editNodes.length; i++) {
                description += '<details>'\
                            + '<summary>'\
                            + editNodes[i].label\
                            + '</summary>'\
                            + '<div style="padding-left: 25">'
                optionsC = {
                    fields: ['id', 'label'],
                    filter: function(item){
                        return item.hidden == false;
                    }
                };
                var connectedNodes = nodes.get(network.getConnectedNodes(editNodes[i].id), optionsC);
                var connectedLabels = [];
                for (var j = 0; j < connectedNodes.length; j++) {
                    connectedLabels.push(connectedNodes[j].label)
                }
                connectedLabels.sort();
                for (var j = 0; j < connectedLabels.length; j++) {
                    description += connectedLabels[j] + '<br>'
                }
                description += '</div></details>'
            }
            if (description) {
                document.getElementById('connections').innerHTML = description
            } else {
                document.getElementById('connections').innerHTML = 'N/A'
            }
            description = '';
            for (var i = 0; i < editNodes.length; i++) {
                description += '<details>'\
                            + '<summary>'\
                            + editNodes[i].label\
                            + '</summary>'\
                            + '<div style="padding-left: 25">'
                var markovNum = editNodes[i].markov
                var markovCluster = nodes.get({
                    filter: function(item) {
                        return (item.markov == markovNum && !item.hidden && item.label != editNodes[i].label);
                    }
                });
                markovCluster = markovCluster.map(a => a.label)
                markovCluster.sort();
                for (var j = 0; j < markovCluster.length; j++) {
                    description += markovCluster[j] + '<br>'
                }
                description += '</div></details>'
            }
        }
        if (description) {
            document.getElementById('markovDesc').innerHTML = description
        } else {
            document.getElementById('markovDesc').innerHTML = 'N/A'
        }
    }
    
    function searchGraph() {
        var searchedValue = document.getElementById('query').value;
        var searchedItems = nodes.get({
            fields: ['id', 'font'],
            filter: function(item) {
                return (item.label.toLowerCase().replace(' ', '_').includes(
                    searchedValue.toLowerCase().replace(' ', '_')
                ));
            }
        });
        let searchedIds = searchedItems.map(a => a.id)
        network.selectNodes(searchedIds);
        if (JSON.stringify(searchedItems) != JSON.stringify(previouslySelected)) {
            for (var i = 0; i < previouslySelected.length; i++) {
                previouslySelected[i].font.background = "transparent";
            }
            nodes.update(previouslySelected);
            for (var i = 0; i < searchedItems.length; i++) {
                searchedItems[i].font.background = "rgb(220,220,220)"
            }
            nodes.update(searchedItems);
            previouslySelected = searchedItems;
        }
        addSelected();
    };
    
    function changeColors() {        
        if (document.getElementById("markov").checked) {
            allNodes = nodes.get({
                fields: ['id', 'color', 'color1']
            })
            for (var i = 0; i < allNodes.length; i++) {
                allNodes[i].color = allNodes[i].color1
            }
        } else if (document.getElementById("clique").checked) {
            allNodes = nodes.get({
                fields: ['id', 'color', 'color0']
            })
            for (var i = 0; i < allNodes.length; i++) {
                allNodes[i].color = allNodes[i].color0
            }
        }
        nodes.update(allNodes);
    };

    function selectGroup() {
        var group = document.getElementById("selectGroup").value;
        if (group == 'selectAll') {
            options = {
                fields: ['id', 'font']
            };
            var groupNodes = nodes.get(options)
        } else if (group == 'selectHidden') {
            options = {
                fields: ['id', 'font'],
                filter: function(item) {
                    return item.hidden
                }
            };
            var groupNodes = nodes.get(options);
        } else if (group == 'selectShown') {
            options = {
                fields: ['id', 'font'],
                filter: function(item) {
                    return !(item.hidden)
                }
            };
            var groupNodes = nodes.get(options)
        } else if (group == 'selectConnected') {
            options = {
                fields: ['id', 'font'],
                filter: function(item){
                    return item.hidden == false;
                }
            };
            var selectedNodes = nodes.get(network.getSelectedNodes(), options);
            var groupNodes = selectedNodes;
            for (var i = 0; i < selectedNodes.length; i++) {
                var connectedNodes = nodes.get(network.getConnectedNodes(selectedNodes[i].id), options);
                groupNodes = [...new Set(groupNodes.concat(connectedNodes))]
            }
        } else if (group == 'selectMarkov') {
            options = {
                fields: ['id', 'font', 'markov']
            };
            var selectedNodes = nodes.get(network.getSelectedNodes(), options);
            var groupNodes = selectedNodes;
            for (var i = 0; i < selectedNodes.length; i++) {
                var markovNum = selectedNodes[i].markov
                var markovCluster = nodes.get({
                    fields: ['id', 'font', 'markov'],
                    filter: function(item) {
                        return (item.markov == markovNum);
                    }
                });
                groupNodes = [...new Set(groupNodes.concat(markovCluster))];
            }
        } else {
            groupNodes = additionalGroups[group];
        };
        groupIds = groupNodes.map(a => a.id);
        network.selectNodes(groupIds);
        addSelected();
        for (var i = 0; i < groupNodes.length; i++) {
            groupNodes[i].font.background = "rgb(220,220,220)"
        }
        nodes.update(groupNodes);
        previouslySelected = [...new Set(previouslySelected.concat(groupNodes))]
    };
    
    function addSelected() {
        var newLine = "\\r\\n";
        var selectedNodes = nodes.get(network.getSelectedNodes());
        var selectedLabels = '';
        for (var i = 0; i < selectedNodes.length; i++) {
            selectedLabels += selectedNodes[i].label + newLine;
        }
        document.getElementById("editNodes").value = selectedLabels.trim();
    }
    
    function hide() {
        var editLabels = document.getElementById("editNodes").value.split("\\n");
        var editIds = [];
        for (var i = 0; i < editLabels.length; i++) {
            editIds.push(labelIdDict[editLabels[i]]);
        }
        editIds = editIds.filter(x => x !== undefined);
        if (editIds.length > 0) {
            var options = {
                fields: ['id', 'hidden', 'physics']
            };
            editNodes = nodes.get(editIds, options);
            for (var i = 0; i < editNodes.length; i++) {
                editNodes[i].hidden = true;
                editNodes[i].physics = false;
            }
            nodes.update(editNodes);
            editEdges = edges.get({
                fields: ['id', 'hidden', 'physics', 'from', 'to'],
                filter: function(item) {
                    return ((editIds.includes(item.from)) || (editIds.includes(item.to)));
                }
            })
            for (var i = 0; i < editEdges.length; i++) {
                editEdges[i].hidden = true;
                editEdges[i].physics = false;
            }
            edges.update(editEdges);
        }
    };
    
    function reveal() {
        var editLabels = document.getElementById("editNodes").value.split("\\n");
        var editIds = [];
        for (var i = 0; i < editLabels.length; i++) {
            editIds.push(labelIdDict[editLabels[i]]);
        }
        editIds = editIds.filter(x => x !== undefined);
        if (editIds.length > 0) {
            var options = {
                fields: ['id', 'hidden', 'physics']
            };
            editNodes = nodes.get(editIds, options);
            for (var i = 0; i < editNodes.length; i++) {
                editNodes[i].hidden = false;
                editNodes[i].physics = true;
            }
            nodes.update(editNodes);
            editEdges = edges.get({
                fields: ['id', 'hidden', 'physics', 'from', 'to'],
                filter: function(item) {
                    return ((editIds.includes(item.from)) || (editIds.includes(item.to)));
                }
            })
            var hiddenNodes = nodes.get({
                filter: function(item) {
                    return (item.hidden)
                }
            });
            var hiddenNodesIds = hiddenNodes.map(a => a.id);
            for (var i = 0; i < editEdges.length; i++) {
                var hide = false;
                if (hiddenNodesIds.includes(editEdges[i].from) || hiddenNodesIds.includes(editEdges[i].to)) {
                    hide = true;
                }
                if (!hide) {
                    editEdges[i].hidden = false;
                    editEdges[i].physics = true;
                }
            }
            edges.update(editEdges);
        }
    };

    function paint() {
        var editLabels = document.getElementById("editNodes").value.split("\\n");
        var editIds = [];
        for (var i = 0; i < editLabels.length; i++) {
            editIds.push(labelIdDict[editLabels[i]]);
        }
        editIds = editIds.filter(x => x !== undefined);
        if (editIds.length > 0) {
            var options = {
                fields: ['id', 'color']
            };
            editNodes = nodes.get(editIds, options);
            var color = document.getElementById("color").value;
            for (var i = 0; i < editNodes.length; i++) {
                editNodes[i].color = color
            }
            nodes.update(editNodes);
        }
    };

    function generateFasta() {
        var editLabels = document.getElementById("editNodes").value.split("\\n");
        var fastaString = '';
        var undefinedFlag = false;
        for (var i = 0; i < editLabels.length; i++) {
            fastaString += '>' + editLabels[i] + '\\n';
            fastaString += fasta[editLabels[i]] + '\\n';
            if (fasta[editLabels[i]] === undefined) {
                undefinedFlag = true;
            }
        }
        if (undefinedFlag) {
            alert("There were undefined fasta sequences! Check your textbox for typos")
        }
        return fastaString;
    }

    function getFasta() {
        var fastaString = generateFasta();
        document.getElementById("modalFasta").style.display = "block";
        document.getElementById("fastaContent").value = fastaString
    };

    function getDesc() {
        describe();
        document.getElementById("modalDesc").style.display = "block";
    };

    function getMSA() {
        document.getElementById("modalMSA").style.display = "block";
    }

    function parseXml(xmlStr) {
        return new window.DOMParser().parseFromString(xmlStr, "text/xml");
    }

    async function fetchClustalId(url, bodyString) {
        const response = await fetch(url, {
            body: bodyString,
            headers: {
                Accept: "text/plain",
                "Content-Type": "application/x-www-form-urlencoded"
            },
            method: "POST"
        })
        return response
    }

    async function getClustalId() {
        var fastaString = generateFasta();
        var emailString = document.getElementById("email").value;
        var bodyString = "email=" \
          + emailString \
          + "&outfmt=fa&sequence=" \
          + fastaString;
        var url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run";
        clustalId = "";
        const response = await fetchClustalId(url, bodyString)
        await response.text().then(r => clustalId = r);
        return response.ok
    }

    function checkStatus() {
        var url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/" \
          + clustalId
        fetch(url, {
          headers: {
              Accept: "text/plain"
          }
        })
            .then(response => response.text().then(r => status = r))
        return status
    }

    function getResults() {
        var url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/"           + clustalId + "/aln-fasta"
        fetch(url, {
            headers: {
                Accept: "text/plain"
            }
        })
            .then(response => response.text().then(r => result = r))
        if (typeof result !== 'undefined') {
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
            element.setAttribute('download', 'result.fa');
            element.style.display = 'none';
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
            return result
        }
    }

    async function runMSA() {
        document.getElementById("clustalStatus").innerHTML = "Please wait...";
        const flag = await getClustalId()
        if (flag) {
            var timeout = 5000;
            var status = "";
            var timerId = setInterval(function(){
                status = checkStatus();
                if (status == "FINISHED") {
                    var result = getResults();
                    if (typeof result !== 'undefined') {
                        var resultObj = result.split('\\n');
                        var msa = {};
                        var currentP = '';
                        for (var i = 0; i < resultObj.length; i++) {
                            if (resultObj[i].charAt(0) == '>') {
                                currentP = resultObj[i]
                            } else {
                                if (!(currentP in msa)) {
                                    msa[currentP] = '';
                                }
                                msa[currentP] += resultObj[i]
                            }
                        }
                        var first = true;
                        for (var key in msa) {
                            if (first) {
                                first = false;
                                var consensus = {};
                                for (var i = 0; i < msa[key].length; i++) {
                                    consensus[i] = {};
                                    consensus[i][msa[key].charAt(i)] = 1;
                                }
                            } else {
                                for (var i = 0; i < msa[key].length; i++) {
                                    currentChar = msa[key].charAt(i);
                                    if (!(currentChar in consensus[i])) {
                                        consensus[i][currentChar] = 0;
                                    }
                                    consensus[i][currentChar] += 1;
                                }
                            }
                        }
                        console.log(consensus);
                        document.getElementById("clustalStatus").innerHTML = ""
                        status = "";
                        clearInterval(timerId);
                    }
                } else if (status == "ERROR" | status == "FAILURE" | status == "NOT_FOUND") {
                    clearInterval(timerId)
                }
                document.getElementById("clustalStatus").innerHTML = "Status: " + status;
            }, timeout)
        } else {
            document.getElementById("clustalStatus").innerHTML = "Status: ERROR";
        }
    }
    
    function closeModal(elemId) {
        document.getElementById(elemId).style.display = "none";
    };

    function assignGroup() {
        var groupName = document.getElementById("groupName").value;
        if (groupName in additionalGroups) {
            alert("This name is already in the drop-down list! Write something else, please");
        } else {
            var editLabels = document.getElementById("editNodes").value.split("\\n");
            var editIds = [];
            for (var i = 0; i < editLabels.length; i++) {
                editIds.push(labelIdDict[editLabels[i]]);
            }
            editIds = editIds.filter(x => x !== undefined);
            options = {
                fields: ['id', 'font'],
            };
            if (editIds.length > 0) {
                var editNodes = nodes.get(editIds, options);
                additionalGroups[groupName] = editNodes;
                var dropDown = document.getElementById("selectGroup");
                var option = document.createElement("option");
                option.text = groupName;
                option.value = groupName;
                dropDown.add(option);
            } else {
                alert("No nodes typed in correctly! Group was not created")
            }
        }
    }

    //nodes.on('*', function (event, properties, senderId) {
    //  console.log('event', event, properties);
    //});

    // BRC code end
''')
            if 'var options, data;' in l:
                f.write('    var options, data;\n')
                f.write('    var additionalGroups = {};\n')
                f.write('    var previouslySelected = [];\n')
                f.write(nodesString.replace('        ', '    '))
                f.write(edgesString.replace('        ', '    '))
                f.write('       ' + fastaJson + ';\n')
            if '</body>' in l:
                f.write('''
<table cellspacing="20">
  <tr>
    <td style="vertical-align:top">
      <h3>Select</h3>
      <div>
        <p>Point and click on desired nodes.<br>Ctrl+click to select multiple nodes.</p>
        <p>
        <input id="query" type="search" name="q" placeholder="Type species/gene names">
        <button id="search" type="button" onclick="return searchGraph()">Select</button>
        </p>
        <em>or</em>
        <p>
        <select id="selectGroup">
          <option value="selectAll">All nodes (hidden & shown)</options>
          <option value="selectHidden">Hidden nodes</options>
          <option value="selectShown">Shown nodes</options>
          <option value="selectConnected">Connected nodes</options>
          <option value="selectMarkov">Markov cluster</options>
        </select>
        <button id="selectGroupButton" type="button" onclick="return selectGroup()">Select group</button>
        </p>
      </div>
    <td style="vertical-align:top">
      <h3>Editing</h3>
      <div>
        <p>
        <h4>Paint nodes by:</h4>
        <input id="clique" name="colors" type="radio" value="clique" checked>
        <label for="clique">Largest maximal cliques</label>
        <br>
        <input id="markov" name="colors" type="radio" value="markov">
        <label for="markov">Markov clustering</label>
        <br>
        <button id="colorChange" onclick="return changeColors()">Change colors</button>
        </p>
      </div>
    </td>
    <td style = "vertical-align:top">
        <h3>Textbox</h4>
      <div>
        <p>
        <textarea id="editNodes" name="editNodes" cols="40" rows="8" placeholder="Lines of node labels" whitespace="pre-wrap"></textarea>
        </p>
        <p>
        <button id="hide" onclick="hide()">Hide</button>
        <button id="reveal" onclick="reveal()">Reveal</button>
        <em>&emsp;or&emsp;</em>
        <input id="color" name="color" placeholder="rgb(r,g,b)" style="float:right">
        <button id="paint" onclick="paint()" style="float:right">Paint</button>
        </p>
        <p>
        <button id="assign" onclick="assignGroup()">Assign group</button>
        <input id="groupName" name="groupName" placeholder="Type your group name" style="float:right">
        </p>
      </div>
    </td>
    <td style="vertical-align:top">
      <h3>Inspect textbox contents</h3>
      <div>
        <p>
        <button id="getFasta" onclick="getFasta()">Get FASTA</button>
        <div id="modalFasta" class="modal">
          <div class="modal-content">
            <span class="close" onclick="closeModal('modalFasta')">&times;</span>
            <h4>FASTA of nodes from the textbox</h4>
            <textarea id="fastaContent" class="fasta"></textarea>
          </div>
        </div>
        </p>
      </div>
      <div>
        <p>
        <button id="getDesc" onclick="getDesc()">Description (visible)</button>
        <div id="modalDesc" class="modal">
          <div class="modal-content">
            <span class="close" onclick="closeModal('modalDesc')">&times;</span>
            <h4>Connections:</h4>
            <div id="connections"></div>
            <h4>Markov Clusters:</h4>
            <div id="markovDesc"></div>
          </div>
        </div>
        </p>
      </div>
      <div>
        <p>
        <button id="getMSA" onclick="getMSA()">Build MSA</button>
        <div id="modalMSA" class="modal">
          <div class="modal-content">
            <span class="close" onclick="closeModal('modalMSA')">&times;</span>
            <p>E-mail:&nbsp;<input id="email" name="email" placeholder="Type your e-mail"></p>
            <button id="clustalButton" type="button" onclick="return runMSA()">Build and download MSA</button>
            <p id="clustalStatus"></p>
            <p id="clustalResults"></p>
          </div>
        </div>
        </p>
      </div>
    </td>
  </tr>
</table>

<script type="text/javascript">
    drawGraph();
    document.getElementById('connections').innerHTML = 'N/A';
    document.getElementById('markovDesc').innerHTML = 'N/A';
</script>\n\n''')
            if (l != nodesString) \
              and (l != edgesString) \
              and (l != commentString) \
              and (not 'var options, data;' in l) \
              and (not 'drawGraph();' in l):
                f.write(l)

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
    with open(os.path.join(rootFolder, 'Results', os.path.splitext(filename)[0] + '.html'),
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
    for filename in os.listdir(inputDir):
        print(filename)
        with open(os.path.join(inputDir, filename), 'r') as oneStrFile:
            mainRefseq = oneStrFile.read().replace('\n', '')
            blast = initialBlast(filename, mainRefseq)
            initBlastList = parseInitialBlast(blast)
            with open(os.path.join(rootFolder, 'Input', filename), 'w') as blastResults:
                blastResults.write('\n'.join(list(dict.fromkeys(initBlastList))))

def runMergeInput():
    '''Run the second step (optional) -
    merge all results into one protein database
    '''
    mergedSet = set()
    for filename in os.listdir(os.path.join(rootFolder, 'Input')):
        with open(os.path.join(rootFolder, 'Input', filename), 'r') as singleFile:
            singleContent = singleFile.read()
        mergedSet = mergedSet | set(singleContent.split('\n'))
        os.remove(os.path.join(rootFolder, 'Input', filename))
    mergedSet.discard('')
    with open(os.path.join(rootFolder, 'Input', 'merged.txt'), 'w') as mergedFile:
        mergedFile.write('\n'.join(mergedSet))

def runMainAnalysis():
    '''Run the third step -
    enrichment of protein database and Blast search
    '''
    for filename in os.listdir(os.path.join(rootFolder, 'Input')):
        proteins = getSequences(filename, dict())
        proteins = getIsoforms(proteins)
        proteins = getSpeciesName(proteins)
        proteins = clearProteins(proteins)
        savePickle(os.path.splitext(filename)[0], proteins, 'Previous_Proteins')
        print(str(datetime.datetime.now()) + ': "proteins" ready')
        blastDict = createBlastDict(proteins, filename)

def runFinalAnalysis():
    '''Run the forth step -
    analysis of Blast results
    '''
    filenameList = os.listdir(inputDir)
    if mergeInput:
        mainFilename = int(input('Choose the main object of the study from the list:\n' + '\n'.join([str(i) + ':' + v for i, v in enumerate(filenameList)]) + '\n'))
        mainFilename = filenameList[mainFilename]
        filenameList = ['merged']
    for filename in filenameList:
        print(filename)
        # proteins need to be refreshed each time we do an analysis
        # else good values are not dropped
        pkl = checkPreviousPickle(filename, 'For_online')
        blastDict = pkl['blastDict']
        pkl = checkPreviousPickle(filename, 'Previous_Proteins')
        proteins = pkl
        transDict, geneDict, greatIso = createDictsForAnalysis(proteins, blastDict)
        if mergeInput:
            filename = mainFilename
        createFastasForTrees(proteins, greatIso, filename)
        with open(os.path.join(inputDir, filename), 'r') as oneStrFile:
            mainRefseq = oneStrFile.read().replace('\n', '').strip()
        for k in proteins.keys():
            if k.split('.')[0] == mainRefseq:
                mainRefseq = k
        mainSpecies = proteins[mainRefseq].species
        mainGene = proteins[mainRefseq].gene 
        graph, maxCliques = createGraph(mainGene, mainSpecies, proteins, geneDict)
        drawGraph(graph, maxCliques, proteins, filename, mainSpecies)
        changeVisJS(filename)
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

if __name__ == '__main__':
    print(str(datetime.datetime.now()) + ': start')
    if preInput:
        runPreInput()
    if mergeInput:
        runMergeInput()
    if doMainAnalysis:
        runMainAnalysis()
    if finalAnalysis:
        runFinalAnalysis()
