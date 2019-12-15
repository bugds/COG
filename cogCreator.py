import sys
import os
import pickle
from io import StringIO
from Bio import SearchIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

rootFolder = sys.path[0]
# 'rootFoler' is a directory that contains:
# /Input        (pairs of input files (good, all) for each protein)
# /Blast_XML    (search results in xml format)
# /Results      (for output)

hitlist_size = 700
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

def checkPreviousProteins(filename):
    for prevName in os.listdir(rootFolder + '/Previous_Proteins'):
        if filename in prevName:
            path = rootFolder + '/Previous_Proteins/' + prevName
            with open(path, 'rb') as f:
                return pickle.load(f)
    return False

def saveProteins(shortName, proteins):
    path = rootFolder + '/Previous_Proteins/' + shortName + '.pkl'
    with open(path, 'wb') as f:
        pickle.dump(proteins, f)


def getSequences(seqFilename, proteins, good=True):
    '''Gets accession numbers from corresponding file

    :param seqFilename: Name of a file with accession numbers
    :param proteins: Dictionary for storing information about proteins
    :param good: Boolean, True for referencial proteins
    :return: Dictionary supplemented with proteins information
    '''
    path = rootFolder + '/Input/' + seqFilename
    seqFile = open(path, 'r')
    line = seqFile.readline()
    while line:
        if not line.replace('\n', '') in proteins:
            proteins[line.replace('\n', '')] = \
                ProteinClass(None, None, line.replace('\n', ''), good)
        line = seqFile.readline()
    return proteins

def getIsoforms(proteins):
    ''' Getting isoforms for all proteins in tree
    
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary supplemented with isoforms
    '''
    record = Entrez.read(Entrez.elink(
        db="gene", 
        dbfrom="protein", 
        id=','.join(proteins.keys())
    ))
    genes = [i['Id'] for i in record[0]['LinkSetDb'][0]['Link']]
    # Safe efetch usage is less than 20 uids at a time (???)
    # IDK why, but if more it sends "An error has occured"
    for i in range(0, 1 + (len(genes) // 20)):
        if i*20 > len(genes):
            efetchAndParse(genes[i*20:len(genes)], proteins)
        else:
            efetchAndParse(genes[i*20:(i+1)*20], proteins) 
    proteins = goodGeneMakesGoodProtein(proteins) 
    return proteins

def efetchAndParse(genesPart, proteins):
    record = Entrez.efetch(
        db="gene", 
        rettype="gene_table", 
        retmode="text", 
        id=','.join(genesPart)
    )

    genesStrings = list()
    line = record.readline()
    if "An error has occured. Please try again later." in line:
        raise TypeError("Something wrong on NCBI side. Try again")
    while '[' in line:
        genesStrings.append(line)
        line = record.readline()
    genesStrings.append(line)

    i=-1
    while line:
        if 'from: ' in line:
            i+=1
            species = genesStrings[i].split('[')[1][:-2]
            gene = genesStrings[i+1].split(',')[0].split(': ')[1]
        if 'Exon table' in line:
            refseq = line.split()[-1]
            if refseq in proteins:
                proteins[refseq].gene = gene
                proteins[refseq].species = species
            elif refseq[1] == 'P':
                proteins[refseq] = ProteinClass(species, gene, refseq, False)
        line = record.readline()

    return proteins

def goodGeneMakesGoodProtein(proteins):
    '''If there are referencial isoforms of some gene,
    all isoforms of this gene must be seen as referencial

    :param goodGenes: List of genes coding referencial isoforms
    :param proteins: Dictionary for storing information about proteins
    :return: "proteins" with all referencial proteins marked
    '''
    goodGenes = [p.gene for p in proteins.values() if p.good]
    for protein in proteins.values():
        if protein.gene in goodGenes:
            protein.good = True
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

def blastSearch(query, species, filename):
    '''Run BLAST, save results of a search to a file and return its contents

    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    '''
    records = NCBIWWW.qblast(
        'blastp',
        'refseq_protein',
        query,
        entrez_query = species,
        hitlist_size = hitlist_size
    )
    xmlPath = rootFolder \
        + '/Blast_XML/' \
        + os.path.splitext(filename)[0] \
        + '.xml'
    xml = open(xmlPath, 'w')
    xml.write(records.getvalue())
    xml.close()
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

def checkBlastDict(filename, blastDict, proteins, iteration):
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
            newBlast = blastSearch(
                '\n'.join(queriesForBlast),
                ' OR '.join(speciesForBlast),
                '{}_iter{}'.format(
                    os.path.splitext(filename)[0], 
                    str(iteration) + '.txt'
                )
            )
        blastDict = createBlastDict(newBlast, blastDict)
        return checkBlastDict(filename, blastDict, proteins, iteration + 1)
    return blastDict

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

def main():
    Entrez.email = email
    for goodFilename in os.listdir(rootFolder + '/Input'):
        if '_good' in goodFilename:
            filename = goodFilename.replace('_good', '')
            shortName = os.path.splitext(filename)[0]
            proteins = checkPreviousProteins(shortName)
            if not proteins:
                proteins = dict()
                proteins = getSequences(goodFilename, proteins)
                proteins = getSequences(filename, proteins, False)
                proteins = getIsoforms(proteins)
                saveProteins(shortName, proteins)
            output = open(rootFolder + '/Results/' + shortName + '.html', 'w')
            blast = checkPreviousBlast(shortName)
            if not blast:
                blast = blastSearch(
                    '\n'.join([p.refseq for p in proteins.values()]),
                    ' OR '.join([p.species for p in proteins.values()]),
                    filename
                )
            blastDict = createBlastDict(blast, dict())
            blastDict = checkBlastDict(filename, blastDict, proteins, 0)
            htmlFull = analyzeBlastDict(blastDict, proteins)
            output.write(htmlFull)
            output.close()

main()