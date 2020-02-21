# Make sure your data is up-to-date;
# Cogger uses NCBI sources which are frequently updated
# Mismatches may cause errors

import sys
import os
import pickle
import subprocess
import time
from io import StringIO
from Bio import SearchIO

rootFolder = sys.path[0]
path2G2R = '/windows1/usr/Boog/gene2refseq/gene2refseq'
path2blastp = '/home/bioinfuser/bioinfuser/ncbi-blast-2.10.0+/bin/blastp'
path2refseq_protein = '/windows1/usr/Boog/BLASTrefseqDB'
path2T2N = '/windows1/usr/Boog/taxid2name/names.dmp'
path2repre = '/windows1/usr/Boog/taxid2name/representatives.txt'

evalueLimit = 0.0001 # Generally 10^-4
qCoverLimit = 0.1 # To get this limit length of a domain of interest must be divided by query length
initBlastTargets = '500'
numThreads = '3'

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
                print(initBlastList[-1])
                print(qCover)
                print(hsp.evalue)
    
    return initBlastList

# This function should be deleted in final version
def checkPreviousPickle(filename, folder):
    '''Takes previously pickled objects or returns "False"

    :param filename: Name of analyzed file
    :param folder: Folder in which pickled objects are contained
    :return: Object or "False" if it does not exist
    '''
    for prevName in os.listdir(rootFolder + folder):
        if filename in prevName:
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
    taxidList = createInputForBlast('.t', taxidList, filename)
    
    time1 = time.time()
    print('start')
    
    bashBlast(
        query=query, 
        out=xmlPath,
        taxidList = taxidList
    )
    
    blast = SearchIO.parse(xmlPath, 'blast-xml')
    
    time2 = time.time()
    print('end (' + str(time2-time1) + ')')
    
    writeInBlastDict(blast, blastDict)
    
    os.remove(query)
    os.remove(species)
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

    subprocess.run(
        [path2blastp, 
        '-db', 'refseq_protein', 
        '-query', query,
        '-outfmt', '5',
        '-out', out,
        '-taxidlist', taxidList,
        '-num_threads', num_threads,
        '-max_target_seqs', max_target_seqs],
        cwd = path2refseq_protein
    )

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
    
def checkBlastDict(p, filename, blastDict, proteins, iteration, previous=set()):
    '''Checks if BLAST found all species in each case

    :param filename: Name of currently explored file
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :param iteration: Number of additional BLAST required currently
    :returns: "blastDict" supported with BLAST results
    '''
    speciesSet = set([p.species for p in proteins.values()])
    speciesForBlast = speciesSet - set(blastDict[p].keys())
    if previous == speciesForBlast:
        for s in speciesForBlast:
            blastDict[p][s] = 'NA'
        return blastDict
    else:
        blastDict = blastSearch(
            p,
            speciesForBlast,
            '{}_iter{}'.format(
                os.path.splitext(filename)[0], 
                str(iteration) + '.nomatter'
            ),
            blastDict
            
        )
        return checkBlastDict(p, filename, blastDict, proteins, iteration + 1,
        speciesForBlast)

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
    htmlString = open(rootFolder + '/htmlStrings.txt', 'r').read()\
        .split('\n')
    htmlString = [line.replace(r'\n', '\n').replace(r'\t', '\t') for line in htmlString]

    htmlPart.write(htmlString[0].format(qSpecies))
    for qGene in qGenes:
        htmlPart.write(htmlString[1].format(
        qGene,
        str(len(gRefseqs) - len(qForward[qGene])),
        len(gRefseqs)
        ))
        for fail in (gRefseqs - qForward[qGene]):
            htmlPart.write(htmlString[2].format(
                fail,
                proteins[fail].species
            ))        
        htmlPart.write(htmlString[3].format(
            str(len(qReverse[qGene]) \
                - len([qR for qR in qReverse[qGene].values() if qR])),
            len(qReverse[qGene])
        ))
        for isoform, success in qReverse[qGene].items():
            htmlPart.write(htmlString[4].format(
                isoform,
                str(len(gSpecies) - len(success)),
                len(gSpecies)
            ))
            for fail in (gSpecies - success):
                htmlPart.write(htmlString[5].format(
                    fail
                ))
            htmlPart.write(htmlString[6])
        htmlPart.write(htmlString[7])
    htmlPart.write(htmlString[8])
    return htmlPart.getvalue()

def mainOffline():
    for filename in os.listdir(rootFolder + '/preInput'):
        print(filename)
        with open(rootFolder + '/preInput/' + filename, 'r') as file:
            blast = initialBlast(filename, file.read().replace('\n', ''))
            parseInitialBlast(blast, qCoverLimit, evalueLimit)
            with open(rootFolder + '/Input/' + os.path.splitext(filename)[0] + '.txt', 'w') as f:
                f.write('\n'.join(list(dict.fromkeys(initBlastList))))


    for filename in os.listdir(rootFolder + '/Input'):
        print(filename)

        proteins = checkPreviousPickle(os.path.splitext(filename)[0], '/Previous_Proteins')

        if not proteins:
            proteins = getSequences(filename, dict())
            start = time.time()            
            proteins = getIsoforms(proteins)
            print(time.time() - start)
            start = time.time()
            proteins = getSpeciesName(proteins)
            print(time.time() - start)
            toDel = list()
            for r in proteins.keys():
                if proteins[r].species == None:
                    toDel.append(r)
                    print('SOMETHING BADD!!!')
                    # print('{}\n{}\n{}\n{}\n{}\n'.format(
                    #     p.species,
                    #     p.taxid,
                    #     p.symbol,
                    #     p.gene,
                    #     p.refseq
                    # ))
            for r in toDel:
                del proteins[r]

            savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')
        
        # print(len(proteins))
        
        # if len(proteins) < 300:
        #     # blast = checkPreviousBlast(os.path.splitext(filename)[0] + '.xml')
        #     # if not blast: ...
        #     blast = blastSearches(
        #         [p.refseq for p in proteins.values()],
        #         [p.taxid for p in proteins.values()],
        #         filename,
        #         blastDict
        #     )
        #     blastDict = createBlastDict(blast, dict())
        #     blastDict = checkBlastDict(filename, blastDict, proteins, 0)
        # else:
        #     blastDict = {}
        #     for p in [p.refseq for p in proteins.values()]:
        #         blastDict = blastSearch(
        #             p, 
        #             [p.taxid for p in proteins.values()], 
        #             filename, 
        #             blastDict
        #         )
        #         blastDict = checkBlastDict(p, filename, blastDict, proteins, 0)

        # savePickle(os.path.splitext(filename)[0], \
        #     {'proteins':proteins, 'blastDict':blastDict}, '/For_online')

mainOffline()