import itertools
import pandas as pd
import json 
import matplotlib.pyplot as plt

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

combinations = list(itertools.product(amino_acids, repeat=3))

pwm_file = '../Benchmarking/IBM_PWM_Seq_1XB_BRD4_0XB.txt'
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

peptides = []
with open('four_mer_Scores.txt', 'r') as infile:
	for line in infile:
		peptides.append(line.strip().split())

old_scores = []
new_scores = []
with open('Seq_1XB_BRD4_0XB_four_mer_Scores.txt', 'a') as outfile:
	outfile.write('Peptide,2XBScore,1XBScore\n')
	for seq,score in peptides:
		old_scores.append(float(score))
		new_score = score_sequence(seq[1:])
		new_scores.append(new_score)
		outfile.write(seq + ',' + score + ',' + str(new_score) + '\n')
		
plt.figure(figsize=(5,4))
plt.scatter(old_scores, new_scores, color='b', edgecolor='black')	
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')

border_thickness = 2
tick_thickness = 2

plt.gca().spines['top'].set_linewidth(border_thickness)
plt.gca().spines['right'].set_linewidth(border_thickness)
plt.gca().spines['bottom'].set_linewidth(border_thickness)
plt.gca().spines['left'].set_linewidth(border_thickness)
plt.gca().xaxis.set_tick_params(width=tick_thickness)
plt.gca().yaxis.set_tick_params(width=tick_thickness)

plt.savefig('Score_4-mers_1XB_2XB.png', dpi=300, bbox_inches='tight')
plt.show()