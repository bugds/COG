import sys
import os
import pickle
import subprocess
import datetime
import networkx
from io import StringIO
from Bio import SearchIO
from copy import deepcopy
from networkx.algorithms.approximation import clique

rootFolder = sys.path[0]
path2G2R = '/windows1/usr/Boog/gene2refseq/gene2refseq' # ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
path2blastp = '/home/bioinfuser/bioinfuser/ncbi-blast-2.10.0+/bin/blastp'
path2refseq_protein = '/windows1/usr/Boog/BLASTrefseqDB'
path2T2N = '/windows1/usr/Boog/taxid2name/names.dmp' # ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
path2repre = '/windows1/usr/Boog/taxid2name/representatives.txt'

evalueLimit = 0.0001 # Generally 10^-4
qCoverLimit = 0.1 # To get this limit length of a domain of interest must be divided by query length
initBlastTargets = '500'
numThreads = '3'
blastChunkSize = 200
orthologyThreshold = 1.0
preInput = False
mergeInput = False
doMainAnalysis = True

# 'rootFolder' is a directory that contains:
# /Input        (pairs of input files (good, all) for each protein)
# /Blast_XML    (search results in xml format)
# /Results      (for output)

class ProteinClass():
    '''Class for proteins
    '''
    def __init__(self, species, taxid, symbol, gene, refseq):
        '''Initialization

        :param species: Species in which proteins are synthetized
        :param gene: Coding gene
        :param refseq: Reference sequence accession number
        :param good: Boolean, if referencial protein - True
        '''
        self.species = species
        self.taxid = taxid
        self.symbol = symbol
        self.gene = gene
        self.refseq = refseq
        self.good = False


def initialBlast(filename, query):
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

def parseInitialBlast(blast, qCoverLimit, evalueLimit):
    initBlastList = list()

    for record in blast:
        queryLen = int(record.seq_len)
        for hit in record:
            for hsp in hit:
                alnSpan = int(hsp.query_span)
                qCover = float("{0:.2f}".format(alnSpan/queryLen))
                if (qCover > qCoverLimit) and (hsp.evalue < evalueLimit):
                    # can't just take hit.accession - does not have accession version
                    substrings = hit.id.split('|')
                    for i in range(len(substrings)):
                        if substrings[i] == 'ref':
                            initBlastList.append(substrings[i+1])
                # print(initBlastList[-1])
                # print(qCover)
                # print(hsp.evalue)
    
    return initBlastList

# This function should be deleted in final version
def checkPreviousPickle(filename, folder):
    '''Takes previously pickled objects or returns "False"

    :param filename: Name of analyzed file
    :param folder: Folder in which pickled objects are contained
    :return: Object or "False" if it does not exist
    '''
    for prevName in os.listdir(rootFolder + folder):
        basePrevName = os.path.splitext(prevName)[0]
        if filename == basePrevName:
            path = rootFolder + folder + '/' + prevName
            with open(path, 'rb') as f:
                return pickle.load(f)
    return False

def savePickle(shortName, toSave, folder):
    '''Saves variables into a pickle file

    :param shortName: Part of the analyzed file name
    :param toSave: Object to save
    :param folder: Folder in which pickled objects are contained
    '''
    path = rootFolder + folder + '/' + shortName + '.pkl'
    with open(path, 'wb') as f:
        pickle.dump(toSave, f)

def getSequences(seqFilename, proteins):
    '''Gets accession numbers from corresponding file

    :param seqFilename: Name of a file with accession numbers
    :param proteins: Dictionary for storing information about proteins
    :param good: Boolean, True for referencial proteins
    :return: Dictionary supplemented with proteins information
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
    ''' Getting isoforms for all proteins in tree
    
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
    return {
            't':header.index('#tax_id'), 
            'g':header.index('GeneID') , 
            'p':header.index('protein_accession.version') , 
            's':header.index('Symbol')
        }

def saveTempSet(toSave, tempSet, proteins, index):
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
            line = f.readline()
    return proteins

def blastSearch(query, speciesList, filename, blastDict):
    '''Run BLAST, save results of a search to a file and return its contents

    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
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
    os.remove(xmlPath)
    
    return blastDict

def createInputForBlast(extension, input, filename):
    with open('{}/Temp/{}{}'.format(rootFolder, filename, extension), 'w') as f:
        if isinstance(input, str):
            f.write(input)
        else:
            f.write('\n'.join(input))
    return '{}/Temp/{}{}'.format(rootFolder, filename, extension)

def bashBlast(query, out, taxidList, db='refseq_protein', outfmt='5', 
    num_threads=numThreads, max_target_seqs='500'):

    blastProcess = subprocess.run(
        [path2blastp, 
        '-db', 'refseq_protein', 
        '-query', query,
        '-outfmt', '5',
        '-out', out,
        '-taxidlist', taxidList,
        '-num_threads', num_threads,
        '-max_target_seqs', max_target_seqs],
        cwd = path2refseq_protein,
        stderr = subprocess.PIPE
    )

    if 'Sequence ID not found' in blastProcess.stderr.decode():
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
            species = hit.description.split('[')[1][:-1]
            if not species in blastDict[record.id]:
                substrings = hit.id.split('|')
                for i in range(len(substrings)):
                    if substrings[i] == 'ref':
                        blastDict[record.id][species] = substrings[i+1]
    return blastDict
    
def checkBlastDict(proteins, filename, blastDict, iteration, previous=[set(), set()]):
    '''Checks if BLAST found all species in each case

    :param filename: Name of currently explored file
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
        blastDict = blastSearch(
            queriesForBlast,
            taxidsForBlast,
            '{}_iter{}'.format(
                os.path.splitext(filename)[0], 
                str(iteration) + '.nomatter'
            ),
            blastDict
        )
        return checkBlastDict(proteins, filename, blastDict, iteration + 1,
        [queriesForBlast, taxidsForBlast])


def goodGeneMakesGoodProtein(proteins, goodGenes):
    for g in goodGenes:
        isoforms = [p.refseq for p in proteins.values() if p.gene == g]
        for i in isoforms:
            proteins[i].good = True
    return proteins

def analyzeBlastDict(blastDict, proteins):
    '''Analysis of a BLAST dictionary

    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: HTML-string containing analysis results
    '''
    htmlFull = ''
    gProteins = [p for p in proteins.values() if p.good]
    gSpecies = set([p.species for p in gProteins])

    for qSpecies in set([p.species for p in proteins.values()]):
        qGenes = set()
        qReverse = dict()
        qForward = dict()

        for qGene in [p.gene for p in proteins.values() if p.species == qSpecies]:
            qGenes.add(qGene)
            qReverse[qGene] = dict() # for qRefseq use qReverse.keys()
            qForward[qGene] = set()
            goodGene = False
            # Reciprocal BLAST
            for qRefseq in [p.refseq for p in proteins.values() if p.gene == qGene]:
                if not proteins[qRefseq].good:
                    qReverse[qGene][qRefseq] = set()
                    if proteins[qRefseq].species in gSpecies:
                        raise ValueError("Same species in both groups - remove one")
                    for s in gSpecies:
                        if blastDict[qRefseq][s] in [p.refseq for p in gProteins]:
                            qReverse[qGene][qRefseq].add(s)
                else:
                    goodGene = True
            # Forward BLAST
            for gRefseq in [p.refseq for p in gProteins]:
                if blastDict[gRefseq][qSpecies] in proteins:
                    if proteins[blastDict[gRefseq][qSpecies]].gene == qGene:
                        qForward[qGene].add(gRefseq)

        if not goodGene:
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
#    htmlString.append('\t<details>\n\t\t<summary>&emsp;Gene id: {}</summary>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{} of {} referencial proteins failed forward BLAST:</summary>\n')
    htmlString.append('\t<details>\n\t\t<summary>&emsp;Gene id: {}</summary>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{}/{} referencial -> this gene. Fails:</summary>\n')
    htmlString.append('\t\t\t\t&emsp;&emsp;&emsp;&emsp;{} [{}]<br>\n')
#    htmlString.append('\t\t</details>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{} of {} isoforms failed to find all referencial proteins in first hit:</summary>\n')
    htmlString.append('\t\t</details>\n\t\t<details>\n\t\t\t<summary>&emsp;&emsp;{}/{} isoforms of this gene -> all referencial isoforms. Fails:</summary>\n')
    htmlString.append('\t\t\t<details>\n\t\t\t\t<summary>&emsp;&emsp;&emsp;{} -> {}/{} referencial isoforms. Fails:</summary>\n')
    htmlString.append('\t\t\t\t\t&emsp;&emsp;&emsp;&emsp;&emsp;{}<br>\n')
    htmlString.append('\t\t\t</details>\n')
    htmlString.append('\t\t</details>\n\t</details>\n')
    htmlString.append('</details>')
    # htmlString = [line.replace(r'\n', '\n').replace(r'\t', '\t') for line in htmlString]

    htmlPart.write(htmlString[0].format(qSpecies))
    for qGene in qGenes:
        htmlPart.write(htmlString[1].format(
        qGene,
#        str(len(gRefseqs) - len(qForward[qGene])),
        len(qForward[qGene]),
        len(gRefseqs)
        ))
        for fail in sorted(list((gRefseqs - qForward[qGene]))):
            htmlPart.write(htmlString[2].format(
                fail,
                proteins[fail].species
            ))
        htmlPart.write(htmlString[3].format(
#            str(len(qReverse[qGene]) \
#                - len([qR for qR in qReverse[qGene].values() if qR])),
            len([qR for qR in qReverse[qGene].values() if qR]),
            len(qReverse[qGene])
        ))
        for isoform, success in qReverse[qGene].items():
            htmlPart.write(htmlString[4].format(
                isoform,
#                str(len(gSpecies) - len(success)),
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

def main():
    if preInput:
        for filename in os.listdir(rootFolder + '/preInput'):
            print(filename)
            with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
                mainRefseq = oneStrFile.read().replace('\n', '')
                blast = initialBlast(filename, mainRefseq)
                initBlastList = parseInitialBlast(blast, qCoverLimit, evalueLimit)
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
            proteins = checkPreviousPickle(
                os.path.splitext(filename)[0], 
                '/Previous_Proteins'
                )

            if not proteins:
                proteins = getSequences(filename, dict())       
                proteins = getIsoforms(proteins)
                proteins = getSpeciesName(proteins)
                toDel = list()
                for r in proteins.keys():
                    if proteins[r].species == None:
                        toDel.append(r)
                        print('SOMETHING BADD!!!')
                for r in toDel:
                    del proteins[r]
            
            savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')
            print(str(datetime.datetime.now()) + ': "proteins" ready')

            blastDict = dict()
            chunksForBlast = dict()
            counter = 0
            for p in proteins.values():
                chunksForBlast[p.refseq] = p
                counter += 1
                if counter >= blastChunkSize:
                    blastDict = blastSearch(
                        [seq.refseq for seq in chunksForBlast.values()],
                        [seq.taxid for seq in proteins.values()],
                        filename,
                        blastDict
                    )
                    print(str(datetime.datetime.now()) + ': Blast search completed')
                    counter = 0
                    chunksForBlast = dict()
                    savePickle('part_' + os.path.splitext(filename)[0], \
                        {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
            blastDict = blastSearch(
                        [seq.refseq for seq in chunksForBlast.values()],
                        [seq.taxid for seq in proteins.values()],
                        filename,
                        blastDict
                    )
            print(str(datetime.datetime.now()) + ': Blast search completed')
            savePickle(os.path.splitext(filename)[0], \
                {'proteins':proteins, 'blastDict':blastDict}, '/For_online')
            
            blastDict = checkBlastDict(proteins, filename, blastDict, 0)
            print(str(datetime.datetime.now()) + ': Blast dictionary checked')
            savePickle('os.path.splitext(filename)[0], \
                {'proteins':proteins, 'blastDict':blastDict}, '/For_online')

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
                            geneDict[g][s] = list(targetGenes.keys())[list(targetGenes.values()).index(max(targetGenes.values()))]

        for filename in os.listdir(rootFolder + '/preInput'):
            print(filename)
            with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
                mainRefseq = oneStrFile.read().replace('\n', '')
            mainSpecies = proteins[mainRefseq].species
            mainGene = proteins[mainRefseq].gene   

            graph = networkx.Graph()
            graph.add_node(mainGene)

            for q in geneDict:
                qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
                if (mainSpecies in geneDict[q]) and (qSpecies in geneDict[mainGene]):
                    if (geneDict[q][mainSpecies] == mainGene) and \
                    (geneDict[mainGene][qSpecies] == q):
                        graph.add_node(q)

            for q in graph.nodes():
                for s in geneDict[q]:
                    for t in graph.nodes():
                        qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
                        tSpecies = [p.species for p in proteins.values() if p.gene == t][0]
                        if (tSpecies in geneDict[q]) and (qSpecies in geneDict[t]):
                            if (q != t) and (geneDict[q][tSpecies] == t) and (geneDict[t][qSpecies] == q):
                                graph.add_edge(q, t)

            maxClique = clique.max_clique(graph)

            # for p in proteins.values():
                # setattr(p, 'good', False)
            proteins = goodGeneMakesGoodProtein(proteins, maxClique)

            refDict = dict()
            for p in proteins.values():
                if p.good:
                    refDict[p.species] = p.gene                
            toDel = set()
            for p in proteins.values():
                if (p.species in refDict.keys()) and (p.gene != refDict[p.species]):
                    toDel.add(p.refseq)
            tempProteins = deepcopy(proteins)
            for refseq in toDel:
                tempProteins.pop(refseq)

            html = analyzeBlastDict(blastDict, tempProteins)
            with open(rootFolder + '/Results/' + os.path.splitext(filename)[0] + '.html', 'w') as out:
                out.write('Original cluster:<br>')
                for gene in maxClique:
                    out.write([p.species for p in proteins.values() if p.gene == gene][0] + ': ' + \
                        [p.symbol for p in proteins.values() if p.gene == gene][0] + ' (' + \
                        gene + ')<br>')
                out.write('<br>Results:<br>')
                out.write(html)

# main()

pkl = checkPreviousPickle('1merged', '/For_online')
proteins = pkl['proteins']
blastDict = pkl['blastDict']

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
                geneDict[g][s] = list(targetGenes.keys())[list(targetGenes.values()).index(max(targetGenes.values()))]

filename = 'CLCNKB.txt'
with open(rootFolder + '/preInput/' + filename, 'r') as oneStrFile:
    mainRefseq = oneStrFile.read().replace('\n', '')
mainSpecies = proteins[mainRefseq].species
mainGene = proteins[mainRefseq].gene   

graph = networkx.Graph()
graph.add_node(mainGene)

for q in geneDict:
    qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
    if (mainSpecies in geneDict[q]) and (qSpecies in geneDict[mainGene]):
        if (geneDict[q][mainSpecies] == mainGene) and \
        (geneDict[mainGene][qSpecies] == q):
            graph.add_node(q)

for q in graph.nodes():
    for s in geneDict[q]:
        for t in graph.nodes():
            qSpecies = [p.species for p in proteins.values() if p.gene == q][0]
            tSpecies = [p.species for p in proteins.values() if p.gene == t][0]
            if (tSpecies in geneDict[q]) and (qSpecies in geneDict[t]):
                if (q != t) and (geneDict[q][tSpecies] == t) and (geneDict[t][qSpecies] == q):
                    graph.add_edge(q, t)

maxClique = clique.max_clique(graph)

proteins = goodGeneMakesGoodProtein(proteins, maxClique)

refDict = dict()
for p in proteins.values():
    if p.good:
        refDict[p.species] = p.gene                
toDel = set()
for p in proteins.values():
    if (p.species in refDict.keys()) and (p.gene != refDict[p.species]):
        toDel.add(p.refseq)
tempProteins = deepcopy(proteins)
for refseq in toDel:
    tempProteins.pop(refseq)

html = analyzeBlastDict(blastDict, tempProteins)
with open(rootFolder + '/Results/' + os.path.splitext(filename)[0] + '.html', 'w') as out:
    out.write('Original cluster:<br>')
    for gene in maxClique:
        out.write([p.species for p in proteins.values() if p.gene == gene][0] + ': ' + \
            [p.symbol for p in proteins.values() if p.gene == gene][0] + ' (' + \
            gene + ')<br>')
    out.write('<br>Results:<br>')
    out.write(html)
