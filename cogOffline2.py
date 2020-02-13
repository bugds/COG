# Make sure your data is up-to-date;
# Cogger uses NCBI sources which are frequently updated
# Mismatches may cause errors

import sys
import os
import pickle
import subprocess
import time
from io import StringIO
from Bio import SearchIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

rootFolder = sys.path[0]
path2G2R = '/windows1/usr/Boog/gene2refseq/g2r.tsv'
path2blastp = '/home/bioinfuser/bioinfuser/ncbi-blast-2.10.0+/bin/blastp'
path2refseq_protein = '/windows1/usr/Boog/BLASTrefseqDB'
path2T2N = '/windows1/usr/Boog/taxid2name/s-names.dmp'

# 'rootFolder' is a directory that contains:
# /Input        (pairs of input files (good, all) for each protein)
# /Blast_XML    (search results in xml format)
# /Results      (for output)

hitlist_size = 300
email = 'bug.dmitrij@gmail.com'

class ProteinClass():
    '''Class for proteins
    '''
    def __init__(self, species, taxid, geneSymbol, gene, refseq):
        '''Initialization

        :param species: Species in which proteins are synthetized
        :param gene: Coding gene
        :param refseq: Reference sequence accession number
        :param good: Boolean, if referencial protein - True
        '''
        self.species = species
        self.taxid = taxid
        self.geneSymbol = geneSymbol
        self.gene = gene
        self.refseq = refseq


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
        line = gene2Refseq.readline().replace('\n', '')

        toSave = False

        while line:
            if tempSet[-1].split('\t')[2] in proteins:
                toSave = True

            if line.split('\t')[1] == tempSet[-1].split('\t')[1]:
                tempSet.append(line)
            else:
                if toSave:
                    for l in tempSet:
                        if l.split('\t')[2] != '-':
                            proteins[l.split('\t')[2]] = ProteinClass(
                                None,
                                l.split('\t')[0], 
                                l.split('\t')[3],
                                l.split('\t')[1], 
                                l.split('\t')[2]
                            )
                    toSave = False
                tempSet = [line]

            line = gene2Refseq.readline().replace('\n', '')

        if toSave:
            for l in tempSet:
                if l.split('\t')[2] != '-':
                    proteins[l.split('\t')[2]] = ProteinClass(
                        None,
                        l.split('\t')[0], 
                        l.split('\t')[3],
                        l.split('\t')[1], 
                        l.split('\t')[2]
                    )
    return proteins

def blastSearches(query, speciesList, filename):
    '''Run BLAST, save results of a search to a file and return its contents

    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    '''
    
    xmlPath = rootFolder \
        + '/Blast_XML/' \
        + os.path.splitext(filename)[0] \
        + '.xml'
    
    with open('{}/Temp/{}{}'.format(rootFolder, filename, '.q'), 'w') as q:
        q.write(query)
    query = '{}/Temp/{}{}'.format(rootFolder, filename, '.q')

    with open('{}/Temp/{}{}'.format(rootFolder, filename, '.s'), 'w') as s:
        s.write('\n'.join(taxidList))
    species = '{}/Temp/{}{}'.format(rootFolder, filename, '.s')
    
    time1 = time.time()
    print('start')
    
    subprocess.run(
        [path2blastp, 
        '-db', 'refseq_protein', 
        '-query', query,
        '-outfmt', '5',
        '-out', '{}/Blast_XML/{}'.format(
            rootFolder, os.path.splitext(filename)[0] + '.xml'
        ),
        '-taxidlist', species,
        '-num_threads', '3'],
        cwd = path2refseq_protein
    )
    
    time2 = time.time()
    print('end (' + str(time2-time1) + ')')
    
    os.remove(query)
    os.remove(species)
    
    return SearchIO.parse(xmlPath, 'blast-xml')

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
    
    with open('{}/Temp/{}{}'.format(rootFolder, filename, '.q'), 'w') as q:
        q.write(query)
    query = '{}/Temp/{}{}'.format(rootFolder, filename, '.q')

    with open('{}/Temp/{}{}'.format(rootFolder, filename, '.s'), 'w') as s:
        s.write('\n'.join(taxidList))
    species = '{}/Temp/{}{}'.format(rootFolder, filename, '.s')
    
    time1 = time.time()
    print('start')
    
    subprocess.run(
        [path2blastp, 
        '-db', 'refseq_protein', 
        '-query', query,
        '-outfmt', '5',
        '-out', '{}/Blast_XML/{}'.format(
            rootFolder, os.path.splitext(filename)[0] + '.xml'
        ),
        '-taxidlist', species,
        '-num_threads', '3'],
        cwd = path2refseq_protein
    )
    blast = SearchIO.parse(xmlPath, 'blast-xml')
    
    time2 = time.time()
    print('end (' + str(time2-time1) + ')')
    
    os.remove(query)
    os.remove(species)
    
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
    
    os.remove(xmlPath)
    
    return blastDict

def createBlastDict(blast, blastDict):
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

def checkBlastDictOld(filename, blastDict, proteins, iteration, previous={'queries':set(),
    'species':set()}):
    '''Checks if BLAST found all species in each case

    :param filename: Name of currently explored file
    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :param iteration: Number of additional BLAST required currently
    :returns: "blastDict" supported with BLAST results
    '''
    speciesSet = set([p.species for p in proteins.values()])
    speciesForBlast = set()
    queriesForBlast = set()
    for record in blastDict.keys():
        lostSpecies = speciesSet - set(blastDict[record].keys())
        if bool(lostSpecies):
            speciesForBlast = speciesForBlast | lostSpecies
            queriesForBlast.add(record)
    if bool(queriesForBlast):
        newBlast = checkPreviousBlast('{}_iter{}'.format(
            os.path.splitext(filename)[0], 
            str(iteration) + '.xml'
        ))
        if not newBlast:
            if (previous['queries'] == queriesForBlast) and \
               (previous['species'] == speciesForBlast):
                for s in speciesForBlast:
                    for q in queriesForBlast:
                        if not s in blastDict[q]:
                            blastDict[q][s] = 'NA'
                return blastDict
            else:
                newBlast = blastSearch(
                    '\n'.join(queriesForBlast),
                    '\n'.join(speciesForBlast),
                    '{}_iter{}'.format(
                        os.path.splitext(filename)[0], 
                        str(iteration) + '.nomatter'
                    )
                )
        blastDict = createBlastDict(newBlast, blastDict)
        return checkBlastDict(filename, blastDict, proteins, iteration + 1,
        {'queries':queriesForBlast, 'species':speciesForBlast})
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

def checkGood(blastDict, proteins):
    ''' Consists of 2 checks:
    oneGenePerSpecies - so there would be no duplication events across
    all good proteins. Prints warning

    trulyGood - checks if BLAST results confirm that all good proteins
    do find each other via reciprocal BLAST. Prints warning

    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: "proteins"
    '''
    goodProteins = [p for p in proteins.values() if p.good]

    # def oneGenePerSpecies(goodProteins=goodProteins):
    speciesSet = set()
    genesSet = set()
    for p in goodProteins:
        if (p.species in speciesSet) and (not p.gene in genesSet):
            print('WARNING!!! Your good set contains paralogs in ' + p.species)
        speciesSet.add(p.species)
        genesSet.add(p.gene)

    # def trulyGood(blastDict=blastDict, proteins=proteins):
    for refseq in proteins:
        if proteins[refseq].good:
            for species in set([p.species for p in proteins.values() if p.good]):
                if blastDict[refseq][species] in proteins:
                    if not proteins[blastDict[refseq][species]].good:
                        print('WARNING!!! Good protein:')
                        print(refseq + ' does not find another good protein in ' \
                            + species)
                        print('Finds ' + blastDict[refseq][species] + ' instead of ' \
                            + str([p.refseq for p in goodProteins if p.species == species]))
                else:
                    print('WARNING!!! Good protein:')
                    print(refseq + ' does not find another good protein in ' \
                        + species)
                    print('Finds ' + blastDict[refseq][species] + ' instead of ' \
                        + str([p.refseq for p in goodProteins if p.species == species]))

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
    Entrez.email = email
    for filename in os.listdir(rootFolder + '/Input'):
        print(filename)

        proteins = checkPreviousPickle(os.path.splitext(filename)[0], '/Previous_Proteins')

        if not proteins:
            proteins = getSequences(filename, dict(), False)
            proteins = getIsoformsOld(proteins)
            savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')
        
        print(len(proteins))
        
        if len(proteins) < 300:
            # blast = checkPreviousBlast(os.path.splitext(filename)[0] + '.xml')
            # if not blast: ...
            blast = blastSearches(
                '\n'.join([p.refseq for p in proteins.values()]),
                [p.taxid for p in proteins.values()],
                filename
            )
            blastDict = createBlastDict(blast, dict())
            blastDict = checkBlastDict(filename, blastDict, proteins, 0)
        else:
            blastDict = {}
            for p in [p.refseq for p in proteins.values()]:
                blastDict = blastSearch(
                    p, 
                    [p.taxid for p in proteins.values()], 
                    filename, 
                    blastDict
                )
                blastDict = checkBlastDict(p, filename, blastDict, proteins, 0)

        savePickle(os.path.splitext(filename)[0], \
            {'proteins':proteins, 'blastDict':blastDict}, '/For_online')

mainOffline()
