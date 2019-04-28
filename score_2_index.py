import os
dirin = 'gene_interaction/score/'
for filename in os.listdir(dirin):
    if filename.endswith(".txt"):
        din = dirin + filename
        i = 0
        idx_score = {}
        for line in open(din):
            if line.strip():
                idx_score[i] = float(line)
                i += 1
        slist = sorted(idx_score, key = lambda key: idx_score[key], reverse = True)
        with open('feat_ranking/' + filename, 'w') as f:
            f.write('\n'.join(map(str,slist)))