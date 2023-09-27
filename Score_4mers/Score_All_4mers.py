import itertools
import pandas as pd
import json 

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

combinations = list(itertools.product(amino_acids, repeat=3))

pwm_file = '../IBM_PWM.txt'
with open(pwm_file, 'r') as file:
	pwm_content = file.read()
	
pwm = json.loads(pwm_content)

## Scoring Function
def score_sequence(sequence):
	score = round(sum([pwm[sequence[i]][i] for i in range(len(sequence))]),3)
	return score

four_mer_scores = []	
for three_mer in combinations:
	three_mer_seq = three_mer[0] + three_mer[1] + three_mer[2]
	four_mer_scores.append(['A' + three_mer_seq, score_sequence(three_mer_seq)])
	four_mer_scores.append(['S' + three_mer_seq, score_sequence(three_mer_seq)])

four_mer_scores = sorted(four_mer_scores, key=lambda x: x[1], reverse=True)

with open('four_mer_Scores.txt', 'a') as outfile:
	for four_mer in four_mer_scores:
		outfile.write(four_mer[0] + ' ' + str(four_mer[1]) + '\n')