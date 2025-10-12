import networkx as nx
import pandas as pd
import numpy as np
class hierarchy:
	G=nx.DiGraph()
	def __init__(self,nodes):
		self.G.add_node('0', depth = 0)
		n = open(nodes,'r')
		for line in n.readlines():
			self.get_nodes(line.strip())

	def get_nodes(self,line):
		node_name=""
		edge_name=""
		nodes = []
		edges = []
		edges.append(['0',line.split('.')[0]])
		for i in line.split('.'):
			node_name+= i
			self.G.add_node(node_name, depth = len(node_name.split('.')))
			node_name+= '.'
		aux = []
		edge_name=line.split('.')[0]
		aux.append(edge_name)
		for i in range(len(line.split('.'))-1):
			edge_name+= '.'
			edge_name+= line.split('.')[i+1]
			aux.append(edge_name)
			edges.append(aux)
			aux = []
			aux.append(edge_name)
		self.G.add_edges_from(edges)
	def stats(self):
		for i in range(self.getHeight()+1):
			print('level' + str(i))
			print(self.getNodesByLevel(i))
			print(len(self.getNodesByLevel(i)))
	def removeNonLeafs(self,df):
		non_leafs = set(self.G.nodes())-set(self.getLeafs())
		df2 = pd.DataFrame.copy(df)
		for i in non_leafs:
			df2 = df2[df2.classification !=  i]
		return(df2)
	def getLeafs(self):
		leafs = []
		for node in self.G.nodes():
			if not self.G.neighbors(node):
				leafs.append(node)
		return(set(leafs))
	def getHeight(self):
		return(max([y['depth'] for x,y in self.G.nodes(data=True)]))
	def getNodesByLevel(self,depth):
		return([x for x,y in self.G.nodes(data=True) if y['depth']==depth])
	def getInnerNodes(self,root,desc):
		for node in self.G.neighbors(root):
			if self.G.neighbors(node):
				desc.append(node)
				self.getInnerNodes(node,desc)
	def getDataFromInnerNodes(self,df):
		desc = []
		self.getInnerNodes('0',desc)
		desc.append('0')
		data = {}
		for node in desc:
			data[node] = self.getDataByParent(node,df)
		return(data)
	def getDataByParent(self,parentNode,df):
		df2 = pd.DataFrame.copy(df)
		if (df2['classification'] == parentNode).any():
			df2.loc[df2['classification']==parentNode, 'classification'] = '#' + parentNode
		for i in self.G.neighbors(parentNode):
			df2.loc[df2.classification.str.startswith(i + "."),'classification'] = str(i)
		for i in set(self.G.nodes())- set(self.G.neighbors(parentNode)):
			df2 = df2[df2.classification!=i]
		df2 = df2.reset_index(drop=True)
		return(df2)
	def getDescendants(self,node):
		return(self.G.neighbors(node))
