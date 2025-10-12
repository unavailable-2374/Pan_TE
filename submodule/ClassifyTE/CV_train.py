import numpy as np
import pandas as pd
import time
import sys
import os
import argparse
from optparse import OptionParser
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import train_test_split


import pickle
from HierStack import hierarchy as hie
from HierStack import model as mo
from HierStack import lcpnb as cl
from HierStack.stackingClassifier import *

def main(dataset_filepath, h, algorithm, cost, gamma, model_filename):
	#------------------------ create file to store overall results --------------------------
	output_filepath = 'CV_Results'
	if not os.path.isdir(output_filepath):
		os.mkdir(output_filepath)
	output_filename = 'CV_Results' +'.txt'
	fd = open(os.path.join(output_filepath, output_filename), 'w')

	model_filepath = 'models'
	if not os.path.isdir(model_filepath):
		os.mkdir(model_filepath)

	pkl_filename = model_filename + str(cost) + "_" + str(gamma) + ".pkl"
	print(f'Training ClassifyTE and saving it to "{pkl_filename}" in "{model_filepath}" directory.')

	# ------------------------- Generate ML Models and test them using 10-Fold CV. -----------------------------
	index = 1
	for index in range(1,11):
		print("Iteration", str(index))
		parent_classifiers = {}
		train_data = pd.read_csv(dataset_filepath + '/' + 'train'+ str(index)+ '.csv', low_memory=False)
		train_data_parent = pd.DataFrame.copy(train_data)
		train_label = pd.DataFrame(train_data['classification'])
		train_data = pd.DataFrame(train_data.drop('classification',axis=1))

		dataInnerNode = h.getDataFromInnerNodes(train_data_parent)

		test_data = pd.read_csv(dataset_filepath + '/' + 'test' +str(index) + '.csv')
		test_label = pd.DataFrame(test_data['classification'])
		test_data = pd.DataFrame(test_data.drop('classification', axis=1))

		m = mo.model(h, algorithm)

		print("----------------Training Started----------------")
		parent_classifiers = m.generate_models(dataInnerNode, index, cost, gamma)

		with open(os.path.join(model_filepath, pkl_filename), "wb") as fm:
			pickle.dump(parent_classifiers, fm)
		fm.close()
		m.test_model(test_data, test_label, parent_classifiers)


		# ---------------- Write fold-wise results to a file ---------------------

		fd.write('FOLD: ' + str(index))
		fd.write('\nPrecision: ' + str(m.precision[-1]))
		fd.write('\nRecall: ' + str(m.recall[-1]))
		fd.write('\nF-measure: ' + str(m.hf[-1]))
		fd.write('\n\nEvaluation per level')

		for i in range(len(m.precision_level[-1])):
			fd.write('\nPrecision Level ' + str(i+1) + ' \t ' + str(m.precision_level[-1][i]))
			fd.write('\nRecall Level ' + str(i+1) + ' \t ' + str(m.recall_level[-1][i]))
			fd.write('\nF-measure Level ' + str(i+1) + ' \t ' + str(m.hf_level[-1][i]))
			fd.write('\n------------------------------------------------------------------------------------------------------')
	    

		fd.write('\n\n================================================================================================================================\n')
	    
		
    	
		fold_output_filename = "fold" + str(index) + "_labels_test_" + ".txt"
		f = open(os.path.join(output_filepath, fold_output_filename) ,'w')
		f.write('true \t predicted\n')
		for i in range(len(m.labels_test)):
			f.write(str(m.labels_test[i][0]) +  '\t ' + str(m.labels_test[i][1]) + '\n')
		f.close()
		index +=1

		# ------------------- Write Average results to a file -------------------
	fd.write('\n\nAVERAGE RESULTS')
	fd.write('\nPrecision :' + str(np.mean(m.precision)) + '\tstd : ' + str(np.std(m.precision)))
	fd.write('\nRecall :' + str(np.mean(m.recall)) + '\tstd : ' + str(np.std(m.recall)))
	fd.write('\nF-measure :' + str(np.mean(m.hf)) + '\tstd : ' + str(np.std(m.hf)))
	fd.write('\n\n')

	for i in range(len(m.precision_level[-1])):
		fd.write('\nPrecision level ' + str(i+1) + ' \t ' + str(np.mean(m.precision_level,axis=0)[i]) + ' \tstd: ' +str(np.std(m.precision_level,axis=0)[i]))
		fd.write('\nRecall level '+ str(i+1) + ' \t ' + str(np.mean(m.recall_level,axis=0)[i]) + ' \tstd: ' + str(np.std(m.recall_level,axis=0)[i]))
		fd.write('\nF-measure level: '+ str(i+1) +  ' \t ' + str(np.mean(m.hf_level,axis=0)[i]) + ' \tstd: ' + str(np.std(m.hf_level,axis=0)[i]))
		fd.write('\n------------------------------------------------------------------------------------------------------')  

	fd.close()


if __name__ == '__main__':

	parser = OptionParser()
	parser.add_option("-n", "--nodes_filepath", dest="node_file", help="Path to node filelist.", default='node.txt')
	parser.add_option("-c", "--c_value", dest="c_value", help="c parameter for SVM")
	parser.add_option("-g", "--gamma_value", dest="gamma_value", help="gamma parameter for SVM")
	parser.add_option("-m", "--model_filename", dest="model_filename", help="output model filename", default='ClassifyTE.pkl')
	
	# Hierarchical classification algorithm can be either:
	# 		non-Leaf Local Classifier per Parent Node (nLCPN)
	# 		Local Classifier per Parent Node and Branch (LCPNB)

	parser.add_option("-a", "--algorithm", dest="algorithm", help="Hierarchical classificationalgorithm LCPNB or nLLCPN.", default='lcpnb')

	(options, args) = parser.parse_args()
	dataset_filepath = "./data/"
	node_filepath = "./nodes/"

	h = hie.hierarchy(node_filepath + options.node_file)
	start_time = time.time()
	main(dataset_filepath, h, options.algorithm, options.c_value, options.gamma_value, options.model_filename)
	total_time = time.time() - start_time
	print("---------------------------------------------Ending Training--------------------------------------------")
	print("\nTotal time elapsed in minutes\t", total_time/60)
