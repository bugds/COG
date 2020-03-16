import os, sys

def readTree(fileName):
    tree = open(fileName, 'r')
    line = tree.readline()
    tree.close()
    line = line.replace('(', '') \
               .replace(')', '') \
               .replace(' ', '_')\
               .replace('\'', '')
    accessionList = line.split(',')
    accessionList = ['_'.join(aN.split('_')[0:2]) for aN in accessionList]
    return accessionList

def main(subtreesFolder, accessionFolder):
    for fileName in os.listdir(subtreesFolder):
        if fileName.endswith('.nwk'):
            accessionNumbers = open(accessionFolder + '/' \
                + '.'.join(fileName.split('.')[:-1]) \
                + '.txt', 'w')
            for aN in readTree(subtreesFolder + '/' + fileName):
                accessionNumbers.write(aN + '\n')
            accessionNumbers.close()

main(sys.path[0] + '/Subtrees_To_Parse',
     sys.path[0] + '/Accession_Numbers')
