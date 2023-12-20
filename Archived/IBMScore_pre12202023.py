import pandas as pd
import json
from io import StringIO
import argparse

def score_sequence(sequence, cutoff=0.5, re_all=False):
	total_score = 0
	num_ibms = 0
	seq_score = {}
	for resi,res in enumerate(sequence):
		if resi < len(sequence) - 3:
			if res == 'A' or res == 'S':
				pep_n = sequence[resi+1:resi+4]
				pep_full = sequence[resi:resi+4]
				score = sum([pwm[pep_n[i]][i] for i in range(len(pep_n))])
				if re_all == True:
					seq_score.update({resi+1:[pep_full,score]})
				if score > cutoff:
					num_ibms += 1
					if re_all == False:
						seq_score.update({resi+1:[pep_full,score]})
					#print(pep_full + '-' + str(resi+1) + ': ' + str(add_score))
				total_score += score
	return total_score, num_ibms, seq_score
	


if __name__ == "__main__":
	# Argument Parser 
	parser = argparse.ArgumentParser(description='Program')
	parser.add_argument('-in', '--Input_FASTA_File', action='store', type=str, required=True,
		help='Name of the text file containing the FASTA sequence of the protein of interest. Carot should not be in same line as sequence, UniProt format preferred.')
	parser.add_argument('-out', '--Output_File', action='store', type=str, required=False,
		help='Name of the text file that will contain the results.')
	parser.add_argument('-cut', '--Cutoff_Value', action='store', type=float, required=False, default=0.5,
		help='Cutoff value above which a peptide is considered to be an IBM.')
	parser.add_argument('-pwm', '--PWSM', action='store', type=str, required=False, default='./IBM_PWM.txt',
		help='Supply a specific position weighted score matrix for scoring.')
	parser.add_argument('-all', '--Return_All', action='store_true', required=False,
		help='Use of this option will return all scored peptides, not just peptides above the cutoff.')

	args = parser.parse_args()

	
	# Import the PWM
	pwm_file = args.PWSM
	with open(pwm_file, 'r') as file:
		pwm_content = file.read()
	
	pwm = json.loads(pwm_content)
		
	## Import FASTA
	input_file = args.Input_FASTA_File 
	if args.Output_File:
		output_file = args.Output_File
	sequence = ''
	protein_id = ''
	skip_next = True
	with open(input_file, 'r') as in_file:
		for line in in_file:
			if skip_next == True:
				if line.startswith('>'):
					protein_id = '>' + line[1:].strip()
					skip_next = False
			elif skip_next == False:
				total_score, num_ibms, seq_scores = score_sequence(line.strip(), cutoff=args.Cutoff_Value, re_all=args.Return_All)
				if args.Output_File:
					with open(output_file, 'a') as out_file:
						out_file.write(protein_id + '\n')
						out_file.write('Total IBM Score: ' + "%.2f" % total_score + '\n')
						out_file.write('Number of IBMs: ' + str(num_ibms) + '\n')

				else:
					print(protein_id)
					print('Total IBM Score: ' + "%.2f" % total_score)
					print('Number of IBMs: ' + str(num_ibms))
						
				for resi in seq_scores.keys():
					if seq_scores[resi][1] > 0.5:
						resi_output = str(resi) + ' ' + seq_scores[resi][0] + ' ' + 'X ' + "%.2f" % seq_scores[resi][1]
					else:
						resi_output = str(resi) + ' ' + seq_scores[resi][0] + ' ' + '- ' + "%.2f" % seq_scores[resi][1]
					if args.Output_File:
						with open(output_file, 'a') as out_file:
							out_file.write(resi_output + '\n')
					else:
						print(resi_output)
				skip_next = True
	

