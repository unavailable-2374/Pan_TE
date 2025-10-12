import networkx as nx
import numpy as np
from HierStack import hierarchy
import pandas as pd
import numpy.linalg
class classification:
	H = None
	def __init__(self,true,h,algorithm):
		self.H = h
		self.algorithm = algorithm
		self.predicted = []
		self.true = []
		self.true.append(str(true))
		self.probs = {}

		while list(self.H.G.predecessors(true)):
			self.true.extend(list(self.H.G.predecessors(true)));
			true = list(self.H.G.predecessors(true))[0]
		self.true.pop(-1)
		self.true.reverse()



	def classify(self, example, parent_classifiers):

	# --------------------- non-Leaf Local Classifier per Parent Node --------------------
		if self.algorithm == "nllcpn":
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
			# self.fix_threshold()

		# ------------- Local CLassifier per Parent Node and Branch	--------------------
		elif self.algorithm == "lcpnb":
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
		# else:
		# 	print("specify classification algorithm.")
		# 	break


	# -------------------------- Return true and predicted values for each test instances ----------------------

	def getPiTi(self):
		return len(set(self.true).intersection(self.predicted))
	def getPi(self):
		return float(len(self.predicted))
	def getTi(self):
		return float(len(self.true))
	def getPiPerLevel(self):
		zi = np.zeros(self.H.getHeight())
		for i in range(1,len(self.predicted)+1):
			zi[i-1] = i
		return zi
	def getTiPerLevel(self):
		ci = np.zeros(self.H.getHeight())
		for i in range(1,len(self.true) + 1):
			ci[i-1] = i
		return ci
	def getPiTiPerLevel(self):
		zici = np.zeros(self.H.getHeight())
		for i in range(1,len(self.true) + 1):
			if len(self.true[:i]) == len(self.predicted[:i]):
				zici[i-1] = len(set(self.true[:i]).intersection(self.predicted[:i]))
		return zici
