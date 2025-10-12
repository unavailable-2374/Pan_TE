import networkx as nx
import numpy as np
from HierStack import hierarchy
from HierStack.stackingClassifier import *
import pandas as pd
import numpy.linalg
class lcpnb:
	H = None
	def __init__(self,h):
		self.H = h
		self.predicted = []
		self.probs = {}


	def classify(self, example, parent_classifiers):
		# ------------- Local CLassifier per Parent Node and Branch	--------------------
		
		branches =  self.H.getNodesByLevel(1)
		predicteds = []
		self.probs.update(dict(zip(parent_classifiers['0'].classes_,parent_classifiers['0'].predict_proba(example.reshape(1,-1))[0])))
		paths = [[branch] for branch in branches]
		while paths:
			path = paths.pop()
			while path[-1] not in self.H.getLeafs():
				if type(parent_classifiers[path[-1]]) is str:
					branch = str(parent_classifiers[path[-1]])
					self.probs[parent_classifiers[path[-1]]] = 1
				else:
					self.probs.update(dict(zip(parent_classifiers[path[-1]].classes_,parent_classifiers[path[-1]].predict_proba(example.reshape(1,-1))[0])))
					non_leafs = set(self.H.G.nodes())-set(self.H.getLeafs())
					non_leafs = set(parent_classifiers[path[-1]].classes_).intersection(non_leafs)
					branch = parent_classifiers[path[-1]].predict(example.reshape(1,-1))[0]
					for not_leafs in non_leafs:
						if not_leafs!= branch:
							new_path = list(path)
							new_path.append(not_leafs)
							paths.append(new_path)
				path.append(branch)
				if branch.startswith('#'):
					break
			predicteds.append(path)
		sum_probs = []
		for o in predicteds:
			sum_probs.append(sum([ self.probs[classification]for classification in o ])/len(o))
		self.predicted = predicteds[np.argmax(sum_probs)]
		
		self.predicted = [ c for c in self.predicted if not c.startswith('#') ]
		return self.predicted
		