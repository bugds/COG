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
# Blastp utility
path2blastp = '/home/bioinfuser/applications/ncbi-blast-2.10.1+/bin/blastp'
# Database with representative taxids
path2refseq_protein = '/home/bioinfuser/data/corgi_files/representative_1'
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
path2T2N = '/home/bioinfuser/data/corgi_files/corgi_oct/names.dmp'
# Representative taxids
path2repre = '/home/bioinfuser/data/corgi_files/corgi_oct/euka_taxids.txt'

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

def initialBlast(filename, query):
    '''Run initial Blast - results will constitute the list of
    candidate-homologs
    :param filename: Name of analyzed file
    :param query: Accession number of a query
    :return: Blast search results in xml-format
    '''
    with open(path2repre, 'r') as repre:
        taxidList = repre.read().split('\n')
    query = createInputForBlast('.q', query, filename)
    xmlPath = rootFolder + '/Blast_XML/' + os.path.splitext(filename)[0] + '.xml'
    taxidList = createInputForBlast('.t', taxidList, filename)
    bashBlast(
        query=query, 
        out=xmlPath,
        taxidList=taxidList,
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

def blastSearch(query, speciesList, filename, blastDict):
    '''Run BLAST, save results of a search to a file and return its contents
    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    :return: Contents of xml-file with Blast results in form of a dictionary
    '''
    xmlPath = rootFolder \
        + '/Blast_XML/' \
        + os.path.splitext(filename)[0] \
        + '.xml'
    query = createInputForBlast('.q', query, filename)
    taxidList = createInputForBlast('.t', speciesList, filename)
    blastNotVoid = bashBlast(
        query=query,
        out=xmlPath,
        taxidList = taxidList
    )
    if blastNotVoid:
        blast = SearchIO.parse(xmlPath, 'blast-xml')
        writeInBlastDict(blast, blastDict)
    os.remove(query)
    os.remove(taxidList)
    if removeXml:
        os.remove(xmlPath)
    return blastDict

def createInputForBlast(extension, inp, filename):
    '''Creating files for Blastp input
    :param extension: ".q" for query, ".t" for taxids
    :param inp: Contents of the input
    :param filename: Name of currently analyzed file
    :return: Path to the temporary file
    '''
    with open('{}/Temp/{}{}'.format(rootFolder, filename, extension), 'w') as f:
        if isinstance(inp, str):
            f.write(inp)
        else:
            f.write('\n'.join(inp))
    return '{}/Temp/{}{}'.format(rootFolder, filename, extension)

def bashBlast(query, out, taxidList, outfmt='5',
  num_threads=numThreads, max_target_seqs='500'):
    '''
    '''
    blastProcess = subprocess.run(
        [path2blastp,
        '-db', 'representative_1',
        '-query', query,
        '-outfmt', '5',
        '-out', out,
#        '-taxidlist', taxidList,
        '-num_threads', num_threads,
        '-max_target_seqs', max_target_seqs],
        stderr = subprocess.PIPE
    )

    if 'Sequence ID not found' in blastProcess.stderr.decode():
        print(blastProcess.stderr.decode())
        exit(1)
        return False
    
    return True

def writeInBlastDict(blast, blastDict):
    '''Create dictionary containing BLAST results
    :param blast: contents of the XML-file with BLAST results
    :return: Dictionary containing BLAST results
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
    
def clearProteins(proteins):
    toDel = list()
    for r in proteins.keys():
        if proteins[r].species == None:
            toDel.append(r)
            print('No info on protein ' + r)
    for r in toDel:
        del proteins[r]
    return proteins

def createBlastDict(proteins, filename):
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
            print(str(datetime.datetime.now()) + ': Blast search completed (chunk ' + str(chunkN) + ')')
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
    print(str(datetime.datetime.now()) + ': Blast search completed (chunk ' + str(chunkN) + ')')
    savePickle(os.path.splitext(filename)[0], \
        {'proteins':proteins, 'blastDict':blastDict}, '/For_online')

    print('Checking Blast dictionary...') 
    blastDict = checkBlastDict(proteins, filename, blastDict, 0)
    print(str(datetime.datetime.now()) + ': Blast dictionary checked')
    savePickle(os.path.splitext(filename)[0], \
        {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
    return blastDict

def checkBlastDict(proteins, filename, blastDict, iteration, previous=[set(), set()]):
    '''Checks if BLAST found all species in each case
    :param filename: Name of currently analyzed file
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :param iteration: Number of additional BLAST required currently
    :returns: "blastDict" supported with BLAST results
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
                print('Running BLAST check with the following parameters:')
                print('Query:')
                print(', '.join([seq.refseq for seq in chunksForBlast.values()]))
                print('Species:')
                print(', '.join(sorted(list(taxidsForBlast))))
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
                print(str(datetime.datetime.now()) + ': Blast search completed (iteration ' + str(iteration) + ', chunk ' + str(chunkN) + ')')
                counter = 0
                chunkN += 1
                chunksForBlast = dict()
                savePickle('part_' + os.path.splitext(filename)[0], \
                    {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
        print('Running BLAST check with the following parameters:')
        print('Query:')
        print(', '.join([seq.refseq for seq in chunksForBlast.values()]))
        print('Species:')
        print(', '.join(sorted(list(taxidsForBlast))))
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
        print(str(datetime.datetime.now()) + ': Blast search completed (iteration ' + str(iteration) + ', chunk ' + str(chunkN) + ')')
        savePickle(os.path.splitext(filename)[0], \
            {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
    
        return checkBlastDict(proteins, filename, blastDict, iteration + 1,
        [queriesForBlast, taxidsForBlast])

def createDictsForAnalysis(proteins, blastDict):
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
    maxLen = clique.node_clique_number(graph, mainGene)
    maxCliques = list()
    for c in clique.find_cliques(graph):
        if len(c) == maxLen:
            maxCliques.append(c)
    return maxCliques


def createGraph(mainGene, mainSpecies, proteins, geneDict):
    graph = networkx.Graph()
    graph.add_node(mainGene)
        
    for q in geneDict:
        qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
        if (mainSpecies in geneDict[q]) and (qSpecies in geneDict[mainGene]):
            if (geneDict[q][mainSpecies] == mainGene) and \
            (geneDict[mainGene][qSpecies] == q):
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

def goodGeneMakesGoodProtein(proteins, goodGenes):
    for g in goodGenes:
        isoforms = [p.refseq for p in proteins.values() if p.gene == g]
        for i in isoforms:
            proteins[i].good = True
    return proteins

def clearProteins2(proteins):
    refDict = dict()
    for p in proteins.values():
        if p.good:
            if p.species in refDict:
                if refDict[p.species] != p.gene:
                    print('Multiple genes are good')
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

def analyzeBlastDict(blastDict, proteins):
    '''Analysis of a BLAST dictionary
    :param blastDict: Dictionary containing BLAST results
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
            #goodGene = False
            # Reciprocal BLAST
            for qRefseq in sorted([p.refseq for p in proteins.values() if p.gene == qGene]):
                #if not proteins[qRefseq].good:
                qReverse[qGene][qRefseq] = set()
                for s in sorted(gSpecies):
                    if blastDict[qRefseq][s] in [p.refseq for p in gProteins]:
                        qReverse[qGene][qRefseq].add(s)
                #else:
                    #goodGene = True
            # Forward BLAST
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

def writeHtml(
    proteins,
    gSpecies,
    gRefseqs,
    qSpecies,
    qGenes,
    qReverse,
    qForward):
    '''Writes BLAST analysis for single gene in form of an HTML-file
    :param proteins: Dictionary for storing information about proteins
    :param gSpecies: Set of good species
    :param gRefseqs: Set of good accession numbers
    :param qSpecies: Species linked to analyzed genes
    :param qGenes: Set of analyzed genes
    :param qReverse: Dictionary of reciprocal BLAST: isoform>species
    :param qForward: Set of species found by forward BLAST
    :return: HTML-string of BLAST analysis for single species
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

def reportHtml(filename, maxClique, proteins, html):
    with open(rootFolder + '/Results/' + os.path.splitext(filename)[0] + '.html', 'w') as out:
        out.write('Original cluster:<br>')
        for gene in sorted(list(maxClique)):
            out.write([p.species for p in proteins.values() if p.gene == gene][0] + ': ' + \
                [p.symbol for p in proteins.values() if p.gene == gene][0] + ' (' + \
                gene + ')<br>')
        out.write('<br>Results:<br>')
        out.write(html)

def main():

    print(str(datetime.datetime.now()) + ': start')
    if preInput:
        for filename in os.listdir(rootFolder + '/preInput'):
            print(filename)
            with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
                mainRefseq = oneStrFile.read().replace('\n', '')
                blast = initialBlast(filename, mainRefseq)
                initBlastList = parseInitialBlast(blast)
                with open(rootFolder + '/Input/' + filename, 'w') as blastResults:
                    blastResults.write('\n'.join(list(dict.fromkeys(initBlastList))))

    if mergeInput:
        mergedSet = set()
        for filename in os.listdir(rootFolder + '/Input'):
            with open(rootFolder + '/Input/' + filename, 'r') as singleFile:
                singleContent = singleFile.read()
            mergedSet = mergedSet | set(singleContent.split('\n'))
        mergedSet.discard('')
        with open(rootFolder + '/Input/merged.txt', 'w') as mergedFile:
            mergedFile.write('\n'.join(mergedSet))

    if doMainAnalysis:
        for filename in os.listdir(rootFolder + '/Input'):
            proteins = getSequences(filename, dict())
            proteins = getIsoforms(proteins)
            proteins = getSpeciesName(proteins)
            proteins = clearProteins(proteins)
            savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')
            print(str(datetime.datetime.now()) + ': "proteins" ready')
            blastDict = createBlastDict(proteins, filename)

    if finalAnalysis:
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
            mainSpecies = proteins[mainRefseq].species
            mainGene = proteins[mainRefseq].gene 
            graph, maxCliques = createGraph(mainGene, mainSpecies, proteins, geneDict)
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

main()
