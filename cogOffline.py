# Make sure your data is up-to-date;
# Cogger uses NCBI sources which are frequently updated
# Mismatches may cause errors

import sys
import os
import pickle
import subprocess
from io import StringIO
from Bio import SearchIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

rootFolder = sys.path[0]
path2G2R = '/windows1/usr/Boog/gene2refseq/g2r.tsv'
path2blastp = '/home/bioinfuser/bioinfuser/ncbi-blast-2.10.0+/bin/blastp'
path2refseq_protein = '/windows1/usr/Boog/BLASTrefseqDB'
path2taxidlist = sys.path[0] + '/taxidlist.txt'

# 'rootFolder' is a directory that contains:
# /Input        (pairs of input files (good, all) for each protein)
# /Blast_XML    (search results in xml format)
# /Results      (for output)

hitlist_size = 300
email = 'bug.dmitrij@gmail.com'

class ProteinClass():
    '''Class for proteins
    '''
    def __init__(self, species, gene, refseq, good):
        '''Initialization

        :param species: Species in which proteins are synthetized
        :param gene: Coding gene
        :param refseq: Reference sequence accession number
        :param good: Boolean, if referencial protein - True
        '''
        self.species = species
        self.gene = gene
        self.refseq = refseq
        self.good = good

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

def getSequences(seqFilename, proteins, good=False):
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
            proteins[line] = ProteinClass(None, None, line, good)
        line = seqFile.readline().replace('\n', '')
    return proteins

def getIsoforms(proteins):
    ''' Getting isoforms for all proteins in tree
    
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with isoforms
    '''
    genes = set()
    
    with open(path2G2R, 'r') as gene2Refseq:
        line = gene2Refseq.readline()
        line = gene2Refseq.readline()
        while line:
            if line.split('\t')[2] in proteins:
                genes.add(line.split('\t')[1])
            line = gene2Refseq.readline()
    
    with open(path2G2R, 'r') as gene2Refseq:
        line = gene2Refseq.readline()
        line = gene2Refseq.readline()
        while line:
            if line.split('\t')[1] in genes:
                proteins[line.split('\t')[2]] = ProteinClass(
                    line.split('\t')[0], 
                    line.split('\t')[1], 
                    line.split('\t')[2], 
                    False
                )
            line = gene2Refseq.readline()
    
    proteins = checkProteins(proteins)
    return proteins

def checkProteins(proteins):
    ''' If outdated proteins are present in the input files,
    NCBI fetches not them, but updated ones. In this situation
    "proteins" contains 2 numbers for the same protein: one old,
    with no gene and species info, but with correct "good" 
    attribute, and new, with no "good" attribute, but with gene
    and species information. This function deletes old and
    supplies new with "good" attribute

    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary with deleted old and supplied new proteins
    '''

    toDel = set()
    toAdd = set()
    for old in proteins.values():
        if old.species == None:
            record = Entrez.read(Entrez.efetch(
                db="protein", 
                rettype='gp',
                retmode='xml',
                id=old.refseq
            ))
            tempor = [r['GBFeature_quals'] for r in record[0]['GBSeq_feature-table'] \
                if r['GBFeature_key'] == 'CDS']
            tempor = [r['GBQualifier_value'] for r in next(iter(tempor)) \
                if r['GBQualifier_name'] == 'db_xref']
            gene = tempor[0].split(':')[1]
            print('Not found: {}'.format(old.refseq))
            toDel.add(old.refseq)
            toAdd.add(gene)
            
    for old in toDel:
        del proteins[old]
    
    with open(path2G2R, 'r') as gene2Refseq:
        line = gene2Refseq.readline()
        line = gene2Refseq.readline()
        while line:
            if line.split('\t')[1] in toAdd:
                proteins[line.split('\t')[2]] = ProteinClass(
                    line.split('\t')[0], 
                    line.split('\t')[1], 
                    line.split('\t')[2], 
                    False
                )
            line = gene2Refseq.readline()
            
    return proteins


def checkPreviousBlast(filename):
    '''Check if there was a previous BLAST with the same name

    :param filename: XML-file containing BLAST information
    :return: False if file was not found or xml contents of a found file
    '''
    for xmlName in os.listdir(rootFolder + '/Blast_XML'):
        if filename in xmlName:
            xmlPath = rootFolder + '/Blast_XML/' + xmlName
            return SearchIO.parse(xmlPath, 'blast-xml')
    return False

def blastSearch(query, species, filename, first=True):
    '''Run BLAST, save results of a search to a file and return its contents

    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    '''
    xmlPath = rootFolder \
        + '/Blast_XML/' \
        + os.path.splitext(filename)[0] \
        + '.xml'
    
    if not first:
        with open('{}/Temp/{}{}'.format(rootFolder, filename, '.q'), 'w') as q:
            q.write(query)
        query = '{}/Temp/{}{}'.format(rootFolder, filename, '.q')
        with open('{}/Temp/{}{}'.format(rootFolder, filename, '.s'), 'w') as s:
            s.write(species)
        species = '{}/Temp/{}{}'.format(rootFolder, filename, '.s')
    else:
        query = '{}/Input/{}'.format(rootFolder, filename)
        species = path2taxidlist
    
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
    
    if not first:
        os.remove(query)
        os.remove(species)
    
    return SearchIO.parse(xmlPath, 'blast-xml')

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

def checkBlastDict(filename, blastDict, proteins, iteration, previous={'queries':set(),
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
                    ),
                    False
                )
        blastDict = createBlastDict(newBlast, blastDict)
        return checkBlastDict(filename, blastDict, proteins, iteration + 1,
        {'queries':queriesForBlast, 'species':speciesForBlast})
    return blastDict

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

        proteins = checkPreviousPickle(os.path.splitext(filename)[0], '/Previous_Proteins')

        if not proteins:
            proteins = getSequences(filename, dict(), False)
            proteins = getIsoforms(proteins)
            savePickle(os.path.splitext(filename)[0], proteins, '/Previous_Proteins')

        blast = checkPreviousBlast(os.path.splitext(filename)[0] + '.xml')
        if not blast:
            blast = blastSearch(
                '\n'.join([p.refseq for p in proteins.values()]),
                '\n'.join([p.species for p in proteins.values()]),
                filename
            )

        blastDict = createBlastDict(blast, dict())
        blastDict = checkBlastDict(filename, blastDict, proteins, 0)

        savePickle(os.path.splitext(filename)[0], \
            {'proteins':proteins, 'blastDict':blastDict}, '/For_online')


mainOffline()
