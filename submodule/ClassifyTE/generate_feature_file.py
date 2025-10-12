#MIT Open License
#Manisha

import pickle
import sys
import os
import subprocess
import shutil
from optparse import OptionParser
from shutil import copy2

def feature_generation(curr_dir1,output_file, feature_dir):
	#copy list.txt to features/ and  features/kanalyze-2.0.0/code
	feature_destpath = os.path.join(curr_dir1, feature_dir)
	#print(f'The feature dest path : {feature_destpath}')
	kanalyzer_destpath = feature_destpath + "/kanalyze-2.0.0/code/"
	kanalyzer_input_destpath = feature_destpath + "/kanalyze-2.0.0/input_data/"
	kanalyzer_output_destpath = feature_destpath + "/kanalyze-2.0.0/output_data/"

	mer2_dir = kanalyzer_output_destpath + '2mer/'
	mer3_dir = kanalyzer_output_destpath + '3mer/'
	mer4_dir = kanalyzer_output_destpath + '4mer/'

	if os.path.isdir(mer2_dir):
		subprocess.run(['rm','-R',kanalyzer_output_destpath + '2mer'])
		os.mkdir(mer2_dir)
	else:
		os.mkdir(mer2_dir)

	if os.path.isdir(mer3_dir):
		subprocess.run(['rm','-R',kanalyzer_output_destpath + '3mer'])
		os.mkdir(mer3_dir)
	else:
		os.mkdir(mer3_dir)

	if os.path.isdir(mer4_dir):
		subprocess.run(['rm','-R',kanalyzer_output_destpath + '4mer'])
		os.mkdir(mer4_dir)
	else:
		os.mkdir(mer4_dir)

	#chmod 775 runKanalyzer_generate_all_features and run this script
	curr_dir = os.getcwd()
	kanalyzer_dir = os.chdir(kanalyzer_destpath)
	curr_dir = os.getcwd() + "/"


	subprocess.run(['chmod', '775', 'runKanalyzer_generate_all_features'])
	subprocess.run(['./runKanalyzer_generate_all_features'])

	change_dir = os.chdir(feature_destpath)
	
	subprocess.run(['javac','KmersFeaturesCollector.java'])
	subprocess.run(['javac','BufferReaderAndWriter.java'])
	subprocess.run(['java', 'KmersFeaturesCollector'])

	if not output_file == "feature_file.csv":
		subprocess.run(['mv', 'feature_file.csv', output_file])
	
	curr_dir2 = os.getcwd()
	change_dir = os.chdir(curr_dir1)
	data_dir = 'data/'
	shutil.copy2(curr_dir2+'/' + output_file, data_dir+output_file)
	change_dir = os.chdir(feature_destpath)
	subprocess.run(['rm', output_file])
	change_dir = os.chdir(curr_dir1)

	

def get_data(fasta_file, feature_dir):
	curr_dir1 = os.getcwd()
	data_filepath = curr_dir1 + "/data/"
	feature_destpath = os.path.join(curr_dir1,feature_dir) 
	kanalyzer_destpath = feature_destpath + "/kanalyze-2.0.0/code/"
	kanalyzer_input_destpath = feature_destpath + "/kanalyze-2.0.0/input_data/"
	kanalyzer_output_destpath = feature_destpath + "/kanalyze-2.0.0/output_data/"

	if os.path.isdir(kanalyzer_input_destpath):
		subprocess.run(['rm','-R',kanalyzer_input_destpath])
		os.mkdir(kanalyzer_input_destpath)
	else:
		os.mkdir(kanalyzer_input_destpath)

	if os.path.isdir(kanalyzer_output_destpath):
		subprocess.run(['rm','-R',kanalyzer_output_destpath])
		os.mkdir(kanalyzer_output_destpath)
	else:
		os.mkdir(kanalyzer_output_destpath)

	sequence = ""
	with open(data_filepath + fasta_file, 'rt') as fp: 
		content = fp.read()
		data = content.split(">")
		#print(data)
		i=0
		for line in data:
			meta = line
			if meta != "":
				i+=1
				of = open(kanalyzer_input_destpath + "seq" + str(i)+".fasta" ,"w")
				of.write(">" + meta)
				of.close()
	fp.close()
	curr_dir2 = os.path.join(os.getcwd(),feature_dir) + "/kanalyze-2.0.0/input_data/"
	change_dir = os.chdir(curr_dir2)
	#files = [f for f in os.listdir('.') if os.path.isfile(f)]
	_, _, files = next(os.walk(curr_dir2))
	files = sorted(files)
	with open("list.txt", "w") as ff:
		for f in files:
			ff.write(f)
			ff.write('\n')
	ff.close()
	shutil.copy2('./list.txt', feature_destpath+"list.txt")
	shutil.copy2('./list.txt', kanalyzer_destpath+"list.txt")
	subprocess.run(['rm', 'list.txt'])


def main():

	
	curr_dir1 = os.getcwd()
	parser = OptionParser()

	parser.add_option("-f", "--filename", dest="filename", help="Name of the fasta file.")
	parser.add_option("-o", "--output", dest="output_filename", help="Name of feature file", default="feature_file.csv")
	parser.add_option("-d", "--featuredir", dest="feature_dir", help="feature directory.", default="features")
	(options, args) = parser.parse_args()
	
	fasta_file = options.filename
	output_file = options.output_filename
	
	feature_dir = os.path.join(os.getcwd(),options.feature_dir) + '/'
	src = os.getcwd() + "/features/"

	if not options.feature_dir == "features":
 		destination = shutil.copytree(src, feature_dir)

	get_data(fasta_file, feature_dir)

	
	#if we want to provide sequence in the command line, uncomment following
	# te_id = sys.argv[1]
	# fasta = sys.argv[2]
	# fasta_filepath = "fasta"
	# if not os.path.isdir(fasta_filepath):
	# 		os.mkdir(fasta_filepath)
	# fasta_filename = str(te_id) + ".fasta"

	# with open(os.path.join(fasta_filepath, fasta_filename), "w") as fd:
	# 	fd.write(">" + str(te_id) + "\n")
	# 	fd.write(fasta)
	# fd.close()

	# with open("list.txt", "w") as f:
	# 	f.write(str(te_id)+ ".fasta")
	# f.close()

	feature_generation(curr_dir1, output_file, feature_dir)
	#subprocess.run(['python', 'evaluate.py','-f','feature_file.csv','-n','node.txt', 'm', ''])
	


if __name__ == '__main__':
	main()