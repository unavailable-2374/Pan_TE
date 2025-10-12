import numpy as np
import pandas as pd
import sys
import os
import csv
import time
import argparse
from optparse import OptionParser
from sklearn.linear_model import LogisticRegression

import pickle
from HierStack import hierarchy as hie
from HierStack import lcpnb as lcpnb
from HierStack import nllcpn as nllcpn
from HierStack.stackingClassifier import *

def getSequenceName(curr_dir, feature_folder):
	input_files_dir = os.path.join(curr_dir + feature_folder) + "/kanalyze-2.0.0/input_data/"
	try:
		_, _, files = next(os.walk(input_files_dir))
	except StopIteration:
		print(f"Error: Directory '{input_files_dir}' does not exist or is empty.")
		print(f"Make sure you have run feature generation for the '{feature_folder.strip('/')}' directory first.")
		sys.exit(1)
	
	fasta_files = [f for f in files if f.endswith('.fasta')]
	fasta_files = sorted(fasta_files)
	seqIDs = []
	
	for file in fasta_files:
		try:
			with open(input_files_dir + file, "r") as f:
				header = f.readline().strip()
				if header.startswith('>'):
					head = header.split(">")
					if len(head) > 1:
						ID = head[1]
						seqIDs.append(ID)
					else:
						print(f"Warning: Invalid header format in file {file}: {header}")
				else:
					print(f"Warning: File {file} does not start with FASTA header")
		except Exception as e:
			print(f"Error reading file {file}: {e}")
			continue
	
	if not seqIDs:
		print(f"Error: No valid FASTA files found in {input_files_dir}")
		sys.exit(1)
		
	return seqIDs

def getLabel(content, predicted):
	if predicted in content:
		label = content[predicted]
		label = label.strip("\n")
		return label


def getCodeLabel(lines):
	content = dict()
	code = []
	label = []	
	for line in lines:
		data = line.split(",")
		code.append(data[0])
		label.append(data[1])
	content = dict(zip(code,label))
	return content

def evaluate_model(test_data, parent_classifiers, algorithm, h):
	labels_evaluate = []
	for i in range(len(test_data)):
		if algorithm == "lcpnb":
			c = lcpnb.lcpnb(h)
		elif algorithm == "nllcpn":
			c = nllcpn.nllcpn(h)
		predicted = c.classify(test_data.iloc[i].values.reshape(1,-1),parent_classifiers)
		labels_evaluate.append(predicted)
	return labels_evaluate


def main(h, data, algorithm, modelname):
	# ------------------------- Generate hierarcical classification for the sequence-----------------------------
	model_filepath = "models/"
	pkl_filename = modelname

	test_data = data.iloc[:, 0:(pow(4,2) + pow(4,3) + pow(4,4))]
	test_label = data.iloc[:,-1]

	parent_classifiers = {}
	#load model
	with open(model_filepath + pkl_filename, 'rb') as fb:
		parent_classifiers = pickle.load(fb)


	print("---------------------------Evaluation Started----------------------------\n")

	labels_test = evaluate_model(test_data, parent_classifiers, algorithm, h )
	
	fb.close()
	return labels_test



if __name__ == '__main__':

	parser = OptionParser()
	parser.add_option("-f", "--filename", dest="filename", help="Name of the feature file.", default="feature_file.csv")
	parser.add_option("-d", "--featuredir", dest="feature_dir", help="feature directory.", default="feature")

	parser.add_option("-n", "--node_file", dest="node_file", help="Path to node filelist.", default="node.txt")
	parser.add_option("-m", "--modelname", dest="modelname", help="Model name")
	
	# Hierarchical classification algorithm can be either:
	# 		non-Leaf Local Classifier per Parent Node (nLCPN)
	# 		Local Classifier per Parent Node and Branch (LCPNB)

	parser.add_option("-a", "--algorithm", dest="algorithm", help="Hierarchical classification algorithm LCPNB or nLLCPN.", default='lcpnb')

	(options, args) = parser.parse_args()

	curr_dir1 = os.getcwd()
	dataset_filepath = curr_dir1 + "/data/"
	node_filepath = curr_dir1 + "/nodes/"
	feature_folder = "/" + options.feature_dir
	seq_names = getSequenceName(curr_dir1, feature_folder)
	

	h = hie.hierarchy(node_filepath + options.node_file)

	start_time = time.time()

	with open(dataset_filepath + options.filename , "r") as csvfile:
		data = pd.read_csv(csvfile, low_memory=False)

	with open("./nodes/tree.txt", "r") as f:
		lines = f.readlines()
		content = getCodeLabel(lines)
	f.close()
	

	hier_label = {}
	hier_label = main(h, data, options.algorithm, options.modelname)

	output_filepath = "output/"
	if not os.path.isdir(output_filepath):
		os.mkdir(output_filepath)
	output_filename = "predicted_out_" + options.feature_dir + ".csv"
	output_txt = "predicted_result_" + options.feature_dir + ".txt"
	f = open(os.path.join(output_filepath, output_filename) ,'w')
	ft = open(os.path.join(output_filepath, output_txt) ,'w')
	f.write("Sequence ID" + "," + "Predicted label" + "\n")
	ft.write("Prediction Results" + "\n")
	
	count = 0
	for k in hier_label:
		name = str(seq_names[count])
		print("Prediction for TE sequence of ID: {}".format(name))
		ft.write("Prediction for TE sequence of ID: {}".format(name))
		j = 1
		predicted_labels = []
		for i in k:
			label = getLabel(content,str(i))	
			predicted_labels.append(str(label))
			print("Predicted level {} : {}".format(str(i), str(label)))
			ft.write("Predicted level {} : {}".format(str(i), str(label)))
			ft.write("\n")
			j = j +1
		ft.write("Final label of TE sequence is {}".format(label))
		ft.write('\n\n')
		ft.write('###############################################################')
		ft.write('\n\n')
		seq_id = name.split(" ")
		print(seq_id[0])
		f.write(seq_id[0])
		f.write(',')
		f.write(label)
		f.write('\n')
		print('\n')
		print('###############################################################')
		print('\n')
		count +=1
	f.close()
	ft.close()
	elapsed_time = time.time() - start_time
	print("\nTotal time elapsed in seconds\t", elapsed_time)
	

