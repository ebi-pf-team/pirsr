import argparse
import csv
import json
import re



def process_row(row, rule):
    sequence_id = row[0]
    model_id = row[1]
    hmm_from = int(row[2])
    hmm_to = int(row[3])
    hmm_align = row[4]
    seq_from = int(row[5])
    seq_to = int(row[6])
    seq_align = row[7]
    global result

    if not sequence_id in result:
        result[sequence_id] = {}

    map = map_hmm_to_seq(hmm_from, hmm_align, seq_align)

    # print(sequence_id + "\t"+model_id)
    rule_sites = []

    for grp in rule['Groups'].keys():

        if model_id in result[sequence_id]:
            next

        pass_count = 0
        
        # print(grp)
        # print(rule['Groups'][grp])

        for pos in rule['Groups'][grp]:
            condition = pos['condition']

            # test = 'ABXXXZ'

            # condition = 'A-B-x(3)-Z'
            condition = re.sub('-', '', condition)
            condition = re.sub('\(', '{', condition)
            condition = re.sub('\)', '}', condition)
            condition = re.sub('x', '.', condition)
            # print('cond:' +condition)



            query_seq = re.sub('-', '', seq_align)
            
            if pos['hmmStart'] < len(map) and pos['hmmEnd'] < len(map):
                target_seq = query_seq[map[pos['hmmStart']]: map[pos['hmmEnd']] + 1]
            else:
                print("Target sequence out of alignment borders for query " +
                      model_id+' on hit '+sequence_id)
                target_seq = ''

            # print("Target: "+target_seq)

            if re.search('\A' + condition + '\Z', target_seq):
                # print("have a match!")
                pass_count += 1
        
        if len(rule['Groups'][grp]) == pass_count:
            # print("PASS!! " + str(pass_count))
            rule_sites.extend(rule['Groups'][grp])

    if rule_sites:
        result[sequence_id][model_id] = {
            'hmmFrom': hmm_from,
            'hmmTo': hmm_to,
            'hmmAlign': hmm_align,
            'seqFrom': seq_from,
            'seqTo': seq_to,
            'seqAlign': seq_align,
            'RuleSites': rule_sites,
            'Scope': rule['Scope'],
        }



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

    ap = argparse.ArgumentParser()

    ap.add_argument("-i", "--query", required=True, help="query tsv input file")
    ap.add_argument("-r", "--rules", required=True, help="processed json rules file")
    ap.add_argument("-o", "--out", required=True, help="output json results file")
    args = vars(ap.parse_args())

    tsv_name = args['query']
    rules_name = args['rules']
    out_file = args['out']
    result = {}


    with open(rules_name) as rulesfile:
        rules_hash = json.load(rulesfile)

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
            
    # print("Global result: ")
    # print(result)

    with open(out_file, 'w') as out_file:
        json.dump(result, out_file, indent=4,)

    print ("Done.")


