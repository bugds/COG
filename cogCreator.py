import sys
import os
from io import StringIO
from Bio import SearchIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

rootFolder = sys.path[0]
# 'rootFoler' is a directory that contains:
# /Fasta        (pairs of fasta files (good, all) for each protein)
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

def getSequences(seqFilename, proteins, good=True):
    '''Gets accession numbers and species from a fasta

    :param seqFilename: Name of a fasta file
    :param proteins: Dictionary for storing information about proteins
    :param good: Boolean, True for referencial proteins
    :return: Dictionary supplemented with proteins information
    '''
    path = rootFolder + '/Fasta/' + seqFilename
    seqFile = open(path, 'r')
    line = seqFile.readline()
    while line:
        if '>' in line:
            refseq = line.split(' ')[0][1:]
            if not (refseq in proteins):
                species = line.split('[')[1][:-2]
                proteins[refseq] = ProteinClass(species, None, refseq, good)
        line = seqFile.readline()
    return proteins

def goodGeneMakesGoodProtein(goodGenes, proteins):
    '''If there are referencial isoforms of some gene,
    all isoforms of this gene must be seen as referencial

    :param goodGenes: List of genes coding referencial isoforms
    :param proteins: Dictionary for storing information about proteins
    :return: "proteins" with all referencial proteins marked
    '''
    for protein in proteins.values():
        if protein.gene in goodGenes:
            protein.good = True
    return proteins

def addIsoform(refseq, i, genes, goodGenes, proteins):
    '''Add isoform of protein that's already in "proteins" dictionary

    :param refseq: Accession number of additional isoform
    :param i: Index of isoform coding gene, listed in "genes"
    :param genes: List of genes that code proteins of interest
    :param goodGenes: List of genes that code referencial proteins
    :param proteins: Dictionary for storing information about proteins
    :return: "proteins" supplemented with single isoform
    '''
    if refseq in proteins:
        proteins[refseq].gene = next(iter(genes[i].values()))
        if proteins[refseq].good == True:
            goodGenes.append(proteins[refseq].gene)
    elif refseq[1] == 'P':
        proteins[refseq] = ProteinClass(
            next(iter(genes[i].keys())),
            next(iter(genes[i].values())),
            refseq,
            False
        )
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
    record = Entrez.efetch(
        db="gene", 
        rettype="gene_table", 
        retmode="text", 
        id=','.join(genes)
    )

    line = record.readline()
    genes = list()
    while ('[' in line):
        species = line.split('[')[1][:-2]
        line = record.readline()
        genes.append({species : line.split(',')[0].split(': ')[1]})
    line = record.readline()

    i = -1
    goodGenes = list()
    while line:
        if 'Reference' in line:
            i+=1
        if 'Exon table' in line:
            proteins = addIsoform(
                line.split()[-1], 
                i, 
                genes, 
                goodGenes, 
                proteins
            )
        line = record.readline()

    proteins = goodGeneMakesGoodProtein(goodGenes, proteins)
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

def blastSearch(query, species, filename, writeToFile=True):
    '''Run BLAST, save results of a search to a file and return its contents

    :param query: String with accession numbers divided by paragraphs
    :param species: String with all species, against which BLAST is performed
    :param filename: Name of original fasta file for saving results of BLAST
    :param writeToFile: Boolean, whether to write to a file or not
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
    return SearchIO.parse(xmlPath, 'blast-xml')

def createBlastDict(blast):
    '''Create dictionary containing BLAST results

    :param blast: contents of the XML-file with BLAST results
    :return: Dictionary containing BLAST results
    '''
    blastDict = {}
    for record in blast:
        blastDict[record.id] = {}
        for hit in record:
            species = hit.description.split('[')[1][:-1]
            if not species in blastDict[record.id]:
                substrings = hit.id.split('|')
                for i in range(len(substrings)):
                    if substrings[i] == 'ref':
                        blastDict[record.id][species] = substrings[i+1]
    return blastDict

def createTemporaryArrays(blastDict, proteins):
    '''Creating temporary arrays for assistance

    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: Dictionary of non-referencial proteins; 
        List of referencial proteins;
        Set of species linked to referencial proteins
    '''
    protDict = dict()
    goodList = list()
    goodSpecies = set()

    for protein in proteins.values():
        if protein.good:
            goodList.append(protein)
            goodSpecies.add(protein.species)
        elif protein.gene in protDict:
            protDict[protein.gene].append(protein)
        else:
            protDict[protein.gene] = [protein]
    return protDict, goodList, goodSpecies

def writeHtml(
    species,
    gene,
    reciprocalPercent,
    forwardPercent,
    failedRecipr,
    forwardIsoforms,
    goodSpecies):
    '''Writes BLAST analysis for single gene in form of an HTML-file

    :param species: Species linked to analyzed gene
    :param gene: Analyzed gene
    :param reciprocalPercent: Percent of isoforms which are top hits 
        of referencial BLAST
    :param forwardPercent: Percent of isoforms, that, when BLASTed
        lead to referencial proteins as top hits
    :failedRecipr: 
    :forwardIsoforms:
    :goodSpecies:
    :return: HTML-string of BLAST analysis for single species
    '''
    htmlPart = StringIO()
    htmlPart.write("""<details>
    <summary>{}</summary>
    <details>
    \t<summary>&emsp;Gene id: {}</summary>
    \t<details>
    \t<summary>&emsp;&emsp;{} {}</summary>""".format(
        species,
        gene,
        'Percent of referencial proteins confirming orthology:',
        reciprocalPercent
    ))
    for failed in failedRecipr:
        htmlPart.write('\n&emsp;&emsp;&emsp;&emsp;{} [{}]<br>'.format(
            failed.refseq,
            failed.species
        ))        
    htmlPart.write("""\t\t\t</details>
    \t\t<details>
    \t\t<summary>&emsp;&emsp;{} {}</summary>""".format(
        'Percent of isoforms finding reference:',
        forwardPercent
    ))
    for isoform, fails in forwardIsoforms.items():
        htmlPart.write('<details>')
        htmlPart.write('<summary>&emsp;&emsp;&emsp;{}: {} {} </summary>'.format(
            isoform,
            '{:.1%}'.format(1 - (len(fails)/len(goodSpecies))),
            'of referencial species confirm orthology'
        ))
        for failed in fails:
            htmlPart.write('\n&emsp;&emsp;&emsp;&emsp;&emsp;{}<br>'.format(
                failed
            ))
        htmlPart.write('</details>')
    htmlPart.write('</details></details></details>')
    return htmlPart.getvalue()

def analyzeBlastDict(blastDict, proteins):
    '''Analysis of a BLAST dictionary

    :param blastDict: Dictionary containing BLAST results
    :param proteins: Dictionary for storing information about proteins
    :return: HTML-string containing analysis results
    '''
    protDict, goodList, goodSpecies = \
        createTemporaryArrays(blastDict, proteins)
    htmlFull = StringIO()

    for gene, isoforms in protDict.items():
        isoformsIds = [i.refseq for i in isoforms]
        isoformsSps = next(iter(isoforms)).species
        failedRecipr = []
        for goodProtein in goodList:
            if not (isoformsSps in blastDict[goodProtein.refseq]):
                failedRecipr.append(goodProtein)
                print('{} not found in BLAST for {}'.format(
                    isoformsSps, goodProtein.refseq
                ))
            elif not (blastDict[goodProtein.refseq][isoformsSps] in isoformsIds):
                failedRecipr.append(goodProtein)
                print('{}\'s first hit is not interest ({})'.format(
                    goodProtein.refseq, isoformsSps
                ))

        forwardIsoforms = {}
        for isoform in isoforms:
            failedGood = []
            for oneGoodSp in goodSpecies:
                if oneGoodSp in blastDict[isoform.refseq]:
                    hit = blastDict[isoform.refseq][oneGoodSp]
                    if hit in proteins:
                        if not proteins[hit].good:
                            failedGood.append(oneGoodSp)
                            print('Unexpected???')
                    else:
                        failedGood.append(oneGoodSp)
                        print('{}\'s first hit not referencial in {}'.format(
                            isoform, oneGoodSp
                        ))
                else:
                    failedGood.append(oneGoodSp)
                    print('{} not found in BLAST for {}'.format(
                        oneGoodSp, isoform
                    ))
            forwardIsoforms[isoform.refseq] = failedGood

        bad = len({k: v for k, v in forwardIsoforms.items() if v != []})

        htmlFull.write(
            writeHtml(
                isoformsSps,
                gene,
                '{:.1%}'.format(1 - (len(failedRecipr)/len(goodList))),
                '{:.1%}'.format(1 - bad/len(isoforms)),
                failedRecipr,
                forwardIsoforms,
                goodSpecies
            )
        )
    return htmlFull

def main():
    Entrez.email = email
    for goodFilename in os.listdir(rootFolder + '/Fasta'):
        if '_good' in goodFilename:
            filename = goodFilename.replace('_good', '')
            shortName = os.path.splitext(filename)[0]
            proteins = dict()
            proteins = getSequences(goodFilename, proteins)
            proteins = getSequences(filename, proteins, False)
            proteins = getIsoforms(proteins)
            output = open(rootFolder + '/Results/' + shortName + '.html', 'w')
            blast = checkPreviousBlast(shortName)
            if not blast:
                blast = blastSearch(
                    '\n'.join([p.refseq for p in proteins.values()]),
                    ' OR '.join([p.species for p in proteins.values()]),
                    filename
                )
            blastDict = createBlastDict(blast)
            htmlFull = analyzeBlastDict(blastDict, proteins)
            output.write(htmlFull.getvalue())
            output.close()

main()