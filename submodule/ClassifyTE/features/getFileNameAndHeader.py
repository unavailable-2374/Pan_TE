import os
import csv

type = {}
data = []

list = "./list.txt"
input_path = "kanalyze-2.0.0/input_data/"
output_path = "./"

with open(list, 'rt') as fp: 
	all_files = fp.readlines()
	
	#print(all_files)
	for file in all_files:
		first_file = file.rstrip("\n")
		#print(first_file)
		with open(input_path + first_file, 'rt') as file:
			header = file.readline()
			print(header)
			head = header.split(">")
			ID = head[1]
			print(ID)		
			type = header.split("|")
			label = type[1].rstrip("\n")
			with open(output_path + "label.csv", 'a') as rep:
				writer = csv.writer(rep, delimiter = ',')
				writer.writerow([first_file,label])
	rep.close()
	file.close()
fp.close()

