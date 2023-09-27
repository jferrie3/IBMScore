import itertools
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import sys

# Create a PWM
def create_pwm(freq_matrix,n):
	pwm = {}
	max_per_pos = [0 for _ in range(n)]
	min_per_pos = [100 for _ in range(n)]
	
	# Define a dictionary to store the background frequency matrix (http://www.nimbios.org/~gross/bioed/webmodules/aminoacid.htm)
	bg_freqs = {'A': 0.074, 'C': 0.033, 'D': 0.059, 'E': 0.058, 'F': 0.040, 'G': 0.074, 'H': 0.029,
		'I': 0.038, 'K': 0.071, 'L': 0.076, 'M': 0.018, 'N': 0.044, 'P': 0.050, 'Q': 0.037,
		'R': 0.042, 'S': 0.081, 'T': 0.062, 'V': 0.068, 'W': 0.013, 'Y': 0.033}
	
	# Compute the PWM
	for aa in freq_matrix:
		pwm[aa] = [0 if freq_matrix[aa][i] == 0 else round((freq_matrix[aa][i] / n) / bg_freqs[aa], 2) for i in range(len(freq_matrix[aa]))]
		for pos,pos_val in enumerate(pwm[aa]):
			if pos_val > max_per_pos[pos]:
				max_per_pos[pos] = pos_val
			if pos_val < min_per_pos[pos]:
				min_per_pos[pos] = pos_val
	max_score = sum(max_per_pos)
	min_score = sum(min_per_pos)
	
	for aa in pwm:
		for pos,pos_val in enumerate(pwm[aa]):
			if pos_val > max_per_pos[pos]:
				max_per_pos[pos] = pos_val
			if pos_val < min_per_pos[pos]:
				min_per_pos[pos] = pos_val
	for aa in pwm:
		pwm.update({aa:[round((x - y) / (n*(z - y)),2) for x, y, z in zip(pwm[aa], min_per_pos, max_per_pos)]})
	return pwm

# Combine Two PWMs
def combine_pwms(pwm1, pwm2, n):
	max_per_pos = [0.0 for _ in range(n)]
	min_per_pos = [1.0 for _ in range(n)]
	
	# Compute the Combined PWM
	for aa in pwm1:
		pwm1.update({aa:[sum(x) for x in zip(pwm1[aa], pwm2[aa])]})
	for aa in pwm1:	
		for pos,pos_val in enumerate(pwm1[aa]):
			if pos_val > max_per_pos[pos]:
				max_per_pos[pos] = pos_val
			if pos_val < min_per_pos[pos]:
				min_per_pos[pos] = pos_val
	max_score = sum(max_per_pos)
	min_score = sum(min_per_pos)


	for aa in pwm1:
		pwm1.update({aa:[round((x - y) / (n*(z - y)),2) for x, y, z in zip(pwm1[aa], min_per_pos, max_per_pos)]})
	return pwm1

# Create a Frequency Matrix
def create_freq_matrix(sequences,n):
	# Define a dictionary to store the frequency matrix
	freq_matrix = {'A': [0 for _ in range(n)], 'C': [0 for _ in range(n)], 'D': [0 for _ in range(n)], 'E': [0 for _ in range(n)], 'F': [0 for _ in range(n)], 'G': [0 for _ in range(n)], 'H': [0 for _ in range(n)],
		'I': [0 for _ in range(n)], 'K': [0 for _ in range(n)], 'L': [0 for _ in range(n)], 'M': [0 for _ in range(n)], 'N': [0 for _ in range(n)], 'P': [0 for _ in range(n)], 'Q': [0 for _ in range(n)],
		'R': [0 for _ in range(n)], 'S': [0 for _ in range(n)], 'T': [0 for _ in range(n)], 'V': [0 for _ in range(n)], 'W': [0 for _ in range(n)], 'Y': [0 for _ in range(n)]}
	
	# Compute the frequency matrix
	for seq in sequences:
		for i, aa in enumerate(seq):
			aa_freq_matrix = freq_matrix[aa]
			aa_freq_matrix[i] +=1
			freq_matrix.update({aa:aa_freq_matrix})
	return freq_matrix

# BLOSUM-ify the PWM
def blosumify(pwm,n):
	## compute the blosumified matrix
	b62 = pd.read_csv('../BLOSUM62.txt', delimiter=r'\s+')
	new_matrix = {'A': [0 for _ in range(n)], 'C': [0 for _ in range(n)], 'D': [0 for _ in range(n)], 'E': [0 for _ in range(n)], 'F': [0 for _ in range(n)], 'G': [0 for _ in range(n)], 'H': [0 for _ in range(n)],
		'I': [0 for _ in range(n)], 'K': [0 for _ in range(n)], 'L': [0 for _ in range(n)], 'M': [0 for _ in range(n)], 'N': [0 for _ in range(n)], 'P': [0 for _ in range(n)], 'Q': [0 for _ in range(n)],
		'R': [0 for _ in range(n)], 'S': [0 for _ in range(n)], 'T': [0 for _ in range(n)], 'V': [0 for _ in range(n)], 'W': [0 for _ in range(n)], 'Y': [0 for _ in range(n)]}
	for parent_aa in pwm:
		blosum_pw_aa = [0 for _ in range (n)]
		for child_aa in pwm:
			blosum_pw_aa = [x + y for x, y in zip(blosum_pw_aa, [i * b62[parent_aa][child_aa] for i in pwm[child_aa]])]
		new_matrix.update({parent_aa:blosum_pw_aa})
	## normalize
	max_per_pos = [0 for _ in range(n)]
	min_per_pos = [0 for _ in range(n)]
	
	for aa in new_matrix:
		for pos,pos_val in enumerate(new_matrix[aa]):
			if pos_val > max_per_pos[pos]:
				max_per_pos[pos] = pos_val
			if pos_val < min_per_pos[pos]:
				min_per_pos[pos] = pos_val

	for aa in new_matrix:
		new_matrix.update({aa:[round((x - y) / (n*(z - y)),2) for x, y, z in zip(new_matrix[aa], min_per_pos, max_per_pos)]})

	return new_matrix

def create_freq_matrix_from_binding(sequences,n):
	# Define a dictionary to store the frequency matrix
	freq_matrix = {'A': [0 for _ in range(n)], 'C': [0 for _ in range(n)], 'D': [0 for _ in range(n)], 'E': [0 for _ in range(n)], 'F': [0 for _ in range(n)], 'G': [0 for _ in range(n)], 'H': [0 for _ in range(n)],
		'I': [0 for _ in range(n)], 'K': [0 for _ in range(n)], 'L': [0 for _ in range(n)], 'M': [0 for _ in range(n)], 'N': [0 for _ in range(n)], 'P': [0 for _ in range(n)], 'Q': [0 for _ in range(n)],
		'R': [0 for _ in range(n)], 'S': [0 for _ in range(n)], 'T': [0 for _ in range(n)], 'V': [0 for _ in range(n)], 'W': [0 for _ in range(n)], 'Y': [0 for _ in range(n)]}
	weight_matrix = {'A': [0 for _ in range(n)], 'C': [0 for _ in range(n)], 'D': [0 for _ in range(n)], 'E': [0 for _ in range(n)], 'F': [0 for _ in range(n)], 'G': [0 for _ in range(n)], 'H': [0 for _ in range(n)],
		'I': [0 for _ in range(n)], 'K': [0 for _ in range(n)], 'L': [0 for _ in range(n)], 'M': [0 for _ in range(n)], 'N': [0 for _ in range(n)], 'P': [0 for _ in range(n)], 'Q': [0 for _ in range(n)],
		'R': [0 for _ in range(n)], 'S': [0 for _ in range(n)], 'T': [0 for _ in range(n)], 'V': [0 for _ in range(n)], 'W': [0 for _ in range(n)], 'Y': [0 for _ in range(n)]}

	# Compute the frequency matrix
	for seq_set in sequences:
		seq = seq_set[0]
		weight = 0.79 - 0.0859*np.log10(seq_set[1])

		for i, aa in enumerate(seq):
			if i > 0:
				aa_freq_matrix = freq_matrix[aa]
				aa_freq_matrix[i-1] += 1
				freq_matrix.update({aa:aa_freq_matrix})
				aa_weight_matrix = weight_matrix[aa]
				aa_weight_matrix[i-1] += weight
				weight_matrix.update({aa:aa_weight_matrix})
	
	for aa in freq_matrix.keys():
		aa_matrix = []
		aa_freq_matrix = freq_matrix[aa]
		aa_weight_matrix = weight_matrix[aa]
		for pos_idx in range(len(aa_freq_matrix)):
			if aa_freq_matrix[pos_idx] > 0:
				aa_matrix.append(aa_weight_matrix[pos_idx] / aa_freq_matrix[pos_idx])
			else:
				aa_matrix.append(0.0)
		freq_matrix.update({aa:aa_matrix})
	
	## normalize
	max_per_pos = [0 for _ in range(n)]
	min_per_pos = [100 for _ in range(n)]
	
	for aa in freq_matrix:
		for pos,pos_val in enumerate(freq_matrix[aa]):
			if pos_val > max_per_pos[pos]:
				max_per_pos[pos] = pos_val
			if pos_val < min_per_pos[pos]:
				min_per_pos[pos] = pos_val

	for aa in freq_matrix:
		freq_matrix.update({aa:[round((x - y) / (n*(z - y)),2) for x, y, z in zip(freq_matrix[aa], min_per_pos, max_per_pos)]})
	
	for aa in freq_matrix:
		aa_set = freq_matrix[aa]
		for pos,score in enumerate(freq_matrix[aa]):
			if score < 0.0:
				aa_set[pos] = 0.0
		freq_matrix.update({aa:aa_set})

	return freq_matrix


## Scoring Function
def score_peptide(peptide, pwm):
	peptide = peptide[1:]
	score = sum([pwm[peptide[i]][i] for i in range(len(peptide))])
	return score
	
	
## Sequences based on enrichment values for XIAP-BIR3 in Eckelman et. al. https://doi.org/10.1038/cdd.2008.6
xiap_sequences = []
all_sequences = []
xiap_bir3_enrichments = [['A', 'I', 'I', 'I', 'K', 'L', 'R', 'R', 'V', 'V', 'V'],
	['A', 'G', 'I', 'K', 'P', 'P', 'P', 'P', 'P', 'P', 'R', 'R', 'S', 'V'],
	['G', 'I', 'I', 'I', 'I', 'K', 'K', 'K', 'P', 'P', 'V']]
xiap_bir3_seq_list = list(itertools.product(*xiap_bir3_enrichments))

for proto_seq in xiap_bir3_seq_list:
	seq = ''
	for aa in proto_seq:
		seq +=aa
	xiap_sequences.append(seq)
	all_sequences.append(seq)
	
	
xiap_freq_mat = create_freq_matrix(xiap_sequences,3)
xiap_pwm = create_pwm(xiap_freq_mat,3)

## Sequences based on enrichment values for cIAP1-BIR3 in Eckelman et. al. https://doi.org/10.1038/cdd.2008.6
ciap_sequences = []

ciap1_bir3_enrichments = [['I', 'I', 'I', 'I', 'K', 'K', 'L', 'Q', 'R', 'R', 'V', 'V', 'V'],
	['A', 'I', 'K', 'L', 'P', 'P', 'P', 'P', 'P', 'R', 'V', 'V'],
	['F', 'F', 'F', 'I', 'I', 'I', 'I', 'I', 'P', 'Y', 'Y', 'Y']]
ciap1_bir3_seq_list = list(itertools.product(*ciap1_bir3_enrichments))

for proto_seq in ciap1_bir3_seq_list:
	seq = ''
	for aa in proto_seq:
		seq +=aa
	ciap_sequences.append(seq)
	all_sequences.append(seq)
	
	
ciap_freq_mat = create_freq_matrix(ciap_sequences,3)
ciap_pwm = create_pwm(ciap_freq_mat,3)


## Sequences based on enrichment value for BIR3-cIAP1 in Kurakin et. al. https://doi.org/10.1002/jmr.809
ser1_sequences = []

ser1_bir3_enrichments = [['R', 'R', 'R', 'R','R', 'R', 'S', 'S', 'V', 'V', 'V', 'M'],
	['V', 'V', 'V', 'V', 'V', 'P', 'P', 'P', 'M', 'T'],
	['W', 'W', 'W', 'W', 'W', 'W', 'W', 'W', 'F']]
ser1_bir3_seq_list = list(itertools.product(*ser1_bir3_enrichments))

for proto_seq in ser1_bir3_seq_list:
	seq = ''
	for aa in proto_seq:
		seq +=aa
	ser1_sequences.append(seq)
	all_sequences.append(seq)

## Generation of PWM from prior literature sequences
all_freq_mat = create_freq_matrix(all_sequences,3)
seq_pwm = create_pwm(all_freq_mat,3)

## BRD4 binding data and generation of BRD4 PWM
brd4_seq_set = [['ATPH', 36.4], ['ATPQ', 148], ['AYAW', 1.16], ['AMAR', 231],
	['ALPQ', 96.8], ['ALPP', 356], ['ATPQ', 201], ['AMAA', 131],
	['AQAQ', 666], ['SNPN', 811], ['SKPP', 301], ['SSPP', 870],
	['STPP', 672], ['SDPY', 356], ['SEPF', 6.10], ['SQPQ', 393], ['AVPI', 2.88]]

brd4_freq = create_freq_matrix_from_binding(brd4_seq_set, 3)
brd4_pwm = create_pwm(brd4_freq,3)

## Construct the final PWM
seq_pwm = blosumify(seq_pwm,3)
#seq_pwm = blosumify(seq_pwm,3)
pwm = seq_pwm
#brd4_pwm = blosumify(brd4_pwm,3)
#pwm = combine_pwms(brd4_pwm, seq_pwm, 3)

## Dump the PWM

json_pwm = json.dumps(pwm)

filename = 'IBM_PWM_Seq_1XB.txt'
with open(filename, 'w') as file:
    file.write(json_pwm)

## Benchmarking
output_array = []
for peptide, binding_affinity in brd4_seq_set:
	score = score_peptide(peptide, pwm)
	output_array.append([peptide, score, binding_affinity])
	
column_names = ['Peptide', 'IBM_Score', 'Binding_Affinity']
df = pd.DataFrame(output_array, columns=column_names)

## Fit and Plot
xvals = np.array(df['Binding_Affinity'].to_list())
yvals = np.array(df['IBM_Score'].to_list())
xvals_log = np.log10(xvals)

slope, intercept, rvalue, pvalue, stderr = linregress(xvals_log, yvals)

x_range = np.linspace(xvals.min(), xvals.max(), 100)
y_fit = slope * np.log10(x_range) + intercept

cutoff_range = np.linspace(0.001, 10000, 100)
y_cutoff = 0 * cutoff_range + 0.5

plt.figure(figsize=(5,4))
plt.scatter(xvals, yvals, color='b', edgecolor='black')
plt.plot(x_range, y_fit, color='black')
plt.plot(cutoff_range, y_cutoff, color='red', linewidth=2, linestyle='--')
plt.title('Sequence + BLOSUM PWM Benchmarking', pad=20, fontsize=14, fontweight='bold')
plt.xlabel('KD (\u03BCM)', fontsize=14, fontweight='bold')
plt.xscale('log')
plt.ylabel('IBM Score', fontsize=14, fontweight='bold')
plt.xlim(0.9, 1100)
plt.ylim(0.35, 0.85)
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


plt.text(1, 0.37, 'y = ' + str(round(slope,4)) + 'log\u2081\u2080(x) + ' + str(round(intercept,2)) + '\nR\u00B2: ' + str(round(rvalue**2,3)), fontsize=12, color='black')

plt.savefig('PWM_Seq_1XB.png', dpi=300, bbox_inches='tight')
plt.show()
