import itertools
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import sys

score_list = []
with open('Seq_1XB_BRD4_0XB_four_mer_Scores.txt') as score_file:
	for line in score_file:
		if 'Peptide' in line:
			continue
		else:
			peptide, score1, score2 = line.strip().split(',')
			score_list.append(float(score2))

score_list = sorted(score_list, reverse=True)

x_range = np.linspace(1, len(score_list), len(score_list))

cutoff_range = np.linspace(-1000, 10*len(score_list), 100)
y_cutoff = 0 * cutoff_range + 0.5

plt.figure(figsize=(5,4))
plt.scatter(x_range, score_list, color='b')
plt.plot(cutoff_range, y_cutoff, color='red', linewidth=2, linestyle='--')
plt.title('IBM Score of 4-mers', pad=20, fontsize=14, fontweight='bold')
plt.xlabel('Peptide Rank', fontsize=14, fontweight='bold')
#plt.xscale('log')
plt.ylabel('IBM Score', fontsize=14, fontweight='bold')
plt.ylim(0.0, 1.0)
plt.xlim(-500, len(score_list)+500)
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

border_thickness = 2
tick_thickness = 2

plt.gca().spines['top'].set_linewidth(border_thickness)
plt.gca().spines['right'].set_linewidth(border_thickness)
plt.gca().spines['bottom'].set_linewidth(border_thickness)
plt.gca().spines['left'].set_linewidth(border_thickness)
plt.gca().xaxis.set_tick_params(width=tick_thickness)
plt.gca().yaxis.set_tick_params(width=tick_thickness)


plt.savefig('Score_4-mers_1XB.png', dpi=300, bbox_inches='tight')
plt.show()
