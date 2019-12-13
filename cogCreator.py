import sys
import os
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
    reciprocalAll,
    forwardFailed,
    forwardAll,
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
    \t<summary>&emsp;&emsp;{} {} {} {}</summary>""".format(
        species,
        gene,
        len(failedRecipr),
        'of',
        reciprocalAll,
        'referencial proteins failed reciprocal BLAST:'
    ))
    for failed in failedRecipr:
        htmlPart.write('\n&emsp;&emsp;&emsp;&emsp;{} [{}]<br>'.format(
            failed.refseq,
            failed.species
        ))        
    htmlPart.write("""\t\t\t</details>
    \t\t<details>
    \t\t<summary>&emsp;&emsp;{} {} {} {}</summary>""".format(
        forwardFailed,
        'of',
        forwardAll,
        'isoforms failed to find all referencial proteins in first hit:'
    ))
    for isoform, fails in forwardIsoforms.items():
        htmlPart.write('<details>')
        htmlPart.write('<summary>&emsp;&emsp;&emsp;{}: {} {} {} {}</summary>'.format(
            isoform,
            len(fails),
            'of',
            len(goodSpecies),
            'referencial proteins not found:'
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
                #print('{} not found in BLAST for {}'.format(
                #    isoformsSps, goodProtein.refseq
                #))
            elif not (blastDict[goodProtein.refseq][isoformsSps] in isoformsIds):
                failedRecipr.append(goodProtein)
                #print('{}\'s first hit is not interest ({})'.format(
                #    goodProtein.refseq, isoformsSps
                #))

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
                        #print('{}\'s first hit not referencial in {}'.format(
                        #    isoform, oneGoodSp
                        #))
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
                len(goodList),
                bad,
                len(isoforms),
                failedRecipr,
                forwardIsoforms,
                goodSpecies
            )
        )
    return htmlFull

def main():
    Entrez.email = email
    for goodFilename in os.listdir(rootFolder + '/Input'):
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
            blastDict = createBlastDict(blast, dict())
            blastDict = checkBlastDict(filename, blastDict, proteins, 0)
            htmlFull = analyzeBlastDict(blastDict, proteins)
            output.write(htmlFull.getvalue())
            output.close()

main()
