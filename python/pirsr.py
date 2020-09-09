import csv
import json


def read_rules_input(file_name):
    # print(file_name)
    with open(file_name) as rulesfile:
        rules_data = json.load(rulesfile)
    # print(rules_data)


    return rules_data

# def read_tsv_input(file_name):
#     print(file_name)
#     with open(file_name) as tsvfile:
#         reader = csv.reader(tsvfile, delimiter='\t')
#         for row in reader:
#             if not bool(row):
#                 break
#             query_id = row[0]
#             print(query_id)
#             process_row(row)


def process_row(row, rule):
    sequence_id = row[0]
    model_id = row[1]
    hmm_from = int(row[2])
    hmm_to = int(row[3])
    hmm_align = row[4]
    seq_from = int(row[5])
    seq_to = int(row[6])
    seq_align = row[7]
    
    map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

    print(sequence_id + "\t"+model_id)

    for grp in rule['Groups'].keys():
        print(grp)
        print(rule['Groups'][grp])




def map_hmm_to_seq(hmm_pos, hmm_seq, seq_seq):
    seq_pos = 0
    map = [0]

    for i in range(0, len(hmm_seq)):
        map[hmm_pos:] = [seq_pos]
        
        if hmm_seq[i:i+1] != '.':
            hmm_pos += 1
        if seq_seq[i:i+1] != '-':
            seq_pos +=1

    return map





if __name__ == '__main__':
    # manual testing example files
    tsv_name = 'align_test.tsv'
    rules_name = 'sr_uru.json'


    rules_hash = read_rules_input(rules_name)

    # print(rules_hash)

    # print(tsv_name)
    with open(tsv_name) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            if not bool(row):
                break
            if row[1] in rules_hash:
                rule = rules_hash[row[1]]
                process_row(row, rule)
            else:
                print('ERROR: nonexistent rule ' + row[1] + ' in rules file')
            



