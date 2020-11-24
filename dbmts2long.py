import csv
import sys
import argparse


def check_args(args=None):
    parser = argparse.ArgumentParser(description="Convert dbMTS to a long format with one line per algorithm/allele/mirna+transcript combination",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help="dbMTS input file.",
                        required=True)
    parser.add_argument('-o', '--output', help="Output file.", required=True)
    arguments = parser.parse_args(args)
    return(arguments)


def convertDBMTS(inputFile, outputFile):
    colDict = {'chr_38': 0, 'pos_38': 1, 'ref_38': 2, 'alt_38': 3,
               'chr': 4, 'pos': 5, 'ref': 6, 'alt': 7,
               'vep_cons': 9, 'vep_trans_id': 10, 'vep_gene_name': 11, 'vep_gene_id': 12, 'vep_cann': 25,
               'dbsnp_150': 31, 'gtex_gene_id.version': 45, 'gtex_tissue': 46,
               'm_ref_score': 131, 'm_ref_mirna': 132, 'm_ref_transcript': 133, 'm_ref_corr_tum': 134, 'm_ref_tum_tiss': 135, 'm_ref_corr_norm': 136, 'm_ref_norm_tiss': 137,
               'm_alt_score': 144, 'm_alt_mirna': 145, 'm_alt_transcript': 146, 'm_alt_corr_tum': 147, 'm_alt_tum_tiss': 148, 'm_alt_corr_norm': 149, 'm_alt_norm_tiss': 150,
               't_ref_score': 161, 't_ref_mirna': 162, 't_ref_transcript': 163, 't_ref_corr_tum': 164, 't_ref_tum_tiss': 165, 't_ref_corr_norm': 166, 't_ref_norm_tiss': 167,
               't_alt_score': 174, 't_alt_mirna': 175, 't_alt_transcript': 176, 't_alt_corr_tum': 177, 't_alt_tum_tiss': 178, 't_alt_corr_norm': 179, 't_alt_norm_tiss': 180,
               'r_ref_score': 191, 'r_ref_mirna': 192, 'r_ref_transcript': 193, 'r_ref_corr_tum': 194, 'r_ref_tum_tiss': 195, 'r_ref_corr_norm': 196, 'r_ref_norm_tiss': 197,
               'r_alt_score': 204, 'r_alt_mirna': 205, 'r_alt_transcript': 206, 'r_alt_corr_tum': 207, 'r_alt_tum_tiss': 208, 'r_alt_corr_norm': 209, 'r_alt_norm_tiss': 210}
    colRetain = list(colDict.values())
    colRetain.sort()
    commonCols = colRetain[0:16]
    newFile = []
    with open(inputFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        newFile.append(['chr_38', 'pos_38', 'ref_38', 'alt_38',
                        'chr', 'pos', 'ref', 'alt',
                        'vep_cons', 'vep_trans_id', 'vep_gene_name', 'vep_gene_id', 'vep_cann',
                        'dbsnp_150', 'gtex_gene_id.version', 'gtex_tissue',
                        'algorithm', 'allele', 'score', 'mirna', 'transcript',
                        'corr_tum', 'tum_tiss', 'corr_norm', 'norm_tiss'])
        for line in reader:
            # discard X chromosome (col 1 & 5) (improperly formatted, and most GWAS don't include it)
            if line[0] == 'X' or line[4] == 'X':
                continue
            lineStart = [line[i] for i in commonCols]  # same for each algorithm & allele
            for alg in ['r', 't', 'm']:
                for allele in ['ref', 'alt']:
                    key = '_'.join((alg, allele))
                    scoreList = line[colDict[key + '_' + 'score']].split(';')
                    mirnaList = line[colDict[key + '_' + 'mirna']].split(';')
                    transcriptList = line[colDict[key + '_' + 'transcript']].split(';')
                    corr_tumList = line[colDict[key + '_' + 'corr_tum']].split('|')
                    tum_tissList = line[colDict[key + '_' + 'tum_tiss']].split('|')
                    corr_normList = line[colDict[key + '_' + 'corr_norm']].split('|')
                    norm_tissList = line[colDict[key + '_' + 'norm_tiss']].split('|')
                    # check that the lists are the same length
                    listOfLists = (scoreList,
                                   mirnaList,
                                   transcriptList,
                                   corr_tumList,
                                   tum_tissList,
                                   corr_normList,
                                   norm_tissList)
                    if len(list(set([len(i) for i in listOfLists]))) != 1:
                        print("WARNING: inconsistent number of entries for score/mirna/transcript/tissues found: \n")
                        print('\t'.join(line))
                        continue
                    # create a new line for each mirna/transcript/score/tissue entry
                    for i, score in enumerate(scoreList):
                        newFile.append(lineStart +
                                       [alg, allele,
                                        score, mirnaList[i], transcriptList[i],
                                        corr_tumList[i], tum_tissList[i],
                                        corr_normList[i], norm_tissList[i]])
    with open(outputFile, 'w') as f:
        for i in newFile:
            f.write('\t'.join(i) + '\n')
        # f.write('\n'.join(['\t'.join(i) for i in newFile]) + '\n')


if __name__ == '__main__':
    cmd_args = check_args(sys.argv[1:])
    convertDBMTS(cmd_args.input, cmd_args.output)
