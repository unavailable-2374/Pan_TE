import networkx as nx
import numpy as np
from HierStack import hierarchy
from HierStack.stackingClassifier import *
import pandas as pd
import numpy.linalg
class nllcpn:
	H = None
	def __init__(self,h):
		self.H = h
		self.predicted = []
		self.probs = {}


	def classify(self, example, parent_classifiers):
		# ------------- Local CLassifier per Parent Node and Branch	--------------------
		branch = parent_classifiers['0'].predict(example.reshape(1,-1))[0]
		self.predicted.append(branch)
		while branch not in self.H.getLeafs():
			if type(parent_classifiers[branch]) is str:
				branch = parent_classifiers[branch]
				self.probs[branch] = 1
			else:
				self.probs[parent_classifiers[branch].predict(example.reshape(1,-1))[0]] = parent_classifiers[branch].predict_proba(example.reshape(1,-1))[0][np.where(parent_classifiers[branch].classes_== parent_classifiers[branch].predict(example.reshape(1,-1))[0])[0]][0]
				branch = parent_classifiers[branch].predict(example.reshape(1,-1))[0]
			if branch.startswith('#'):
				break
			self.predicted.append(branch)
			
		return self.predicted
		