import pandas
import subprocess

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

gene2refseq_path = '/home/bioinfuser/data/corgi_files/corgi_2021/gene2refseq'
blastdbcmd_path = '/home/bioinfuser/applications/ncbi-blast-2.10.1+/bin/blastdbcmd'
newDatabase_batch = '/home/bioinfuser/data/corgi_files/diss_batch'
newDatabase_fasta = '/home/bioinfuser/data/corgi_files/diss.fasta'
newDatabase_name = '/home/bioinfuser/data/corgi_files/diss'
refseq_path = '/home/bioinfuser/data/refseq_protein_2021'

email = "bug.dmitrij@gmail.com"

with open('/home/bioinfuser/data/corgi_files/corgi_2021/diss_taxids.txt', 'r') as inpObj:
    taxids = [t.replace('\n', '') for t in inpObj.readlines()]

geneMapDF = pandas.read_csv(gene2refseq_path, sep = '\t', usecols=['#tax_id', 'GeneID', 'protein_accession.version'])
geneMapDF = geneMapDF[geneMapDF['#tax_id'].isin(taxids)]
geneMapDF.drop_duplicates(inplace = True)

with open(newDatabase_batch, 'w') as inpObj:
    inpObj.write('\n'.join([refseq for refseq in geneMapDF['protein_accession.version'] if refseq != '-']))

fastaSearch = subprocess.run(
    [
        blastdbcmd_path,
        '-entry_batch', newDatabase_batch,
        '-db', 'refseq_protein',
        '-out', newDatabase_fasta
    ],
    cwd = refseq_path,
    stderr = subprocess.PIPE
)

errRefseqs_all = fastaSearch.stderr.decode().split('\n')


errRefseqs_all = [i.split(' ')[-1] for i in errRefseqs_all if "Skipped" in i]

with open('/home/bioinfuser/data/corgi_files/errors.txt', 'w') as o:
    i = 0
    for l in errRefseqs_all:
        o.write(l + '\n')
        i += 1
    print('Sequences with errors: ' + str(i))

from Bio import Entrez
Entrez.email = email

errRefseqs_list = list()
counter = 0
step = 600
while counter < len(errRefseqs_all):
    errRefseqs_list.append(errRefseqs_all[counter : counter+step])
    counter += step
errRefseqs_list.append(errRefseqs_all[counter:])

handle_list = list()
for errRefseqs in errRefseqs_list:
    handle_list.append(Entrez.efetch(db="protein", id = ','.join(errRefseqs), rettype="fasta", retmode="text"))

record_list = list()
for handle in handle_list:
    record_list.append(handle.read())

with open(newDatabase_fasta, 'a') as inpObj:
    for record in record_list:
        if not 'Supplied id parameter is empty.' in record:
            inpObj.write(record)
        else:
            print(record)

#   The synonymical fasta sequences are split on several independent sequences,
# each with its own id. In this way blast database will be parsed for queries
# correctly.

from Bio import SeqIO

refseq_set = set()
val = set(geneMapDF['protein_accession.version'].values)

with open(newDatabase_fasta + '_nodup.fasta', 'w') as inpObj:
    newDatabase_fasta_object = SeqIO.parse(newDatabase_fasta, "fasta")
    for record in newDatabase_fasta_object:
        synonyms = [s for s in record.description.split('>') if s != '']
        toDel = set()
        for i in range(0, len(synonyms)):
            if not ((synonyms[i][1:3] == 'P_') and (isfloat('1' + synonyms[i].split(' ')[0][3:]))):
                synonyms[i-1] = '>'.join([synonyms[i-1], synonyms[i]])
                toDel.add(i)
        j = 0
        for i in toDel:
            synonyms.pop(i-j)
            j += 1

        for s in synonyms:
            if not s in refseq_set:
                refseq_set.add(s)
                synonym_record = record
                synonym_record.description = s
                synonym_record.id = s.split(' ')[0]
                synonym_record.name = synonym_record.id
                if synonym_record.id in val:
                    SeqIO.write(synonym_record, inpObj, "fasta")

makeblastdb_path = '/home/bioinfuser/applications/ncbi-blast-2.10.1+/bin/makeblastdb'

makeDB = subprocess.run(
    [
        makeblastdb_path,
        '-dbtype', 'prot',
        '-in', newDatabase_fasta + '_nodup.fasta',
        '-title', newDatabase_name,
        '-parse_seqids',
        '-out', newDatabase_name,
        '-max_file_sz', '4GB'
    ],
    stderr = subprocess.PIPE
)
