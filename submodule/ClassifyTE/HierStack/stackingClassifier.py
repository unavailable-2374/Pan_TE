from sklearn.base import *
import numpy as np
import pandas as pd
import os
import pickle
from sklearn.model_selection import cross_val_predict,LeaveOneOut

from copy import deepcopy


def get_probabilities(model, X):
	pred = model.predict_proba(X)	
	return pred

def get_predictions(model, X):
	pred = model.predict(X)	
	return pred

def get_classes(model, X):
	classes = model.classes_
	return classes


class StackingClassifier:
	base_models = []
	meta_classifier = []
	copy_base_models = []
	copy_meta_model = None


	def __init__(self, classifiers = [], meta_classifier = []):
		self.base_models = classifiers
		self.meta_classifier = meta_classifier
		self.classes_ = []
		self.copy_base_models = []
		self.meta_classifier

	
	def fit_base_models(self, X, y):
		cv_probabilities = []
		index = 1 
		for i, clf in enumerate(self.base_models):
			print('INSIDE 50FOLD CV')
			print(f'---------------------------fitting base model {index}-------------------------')
			clf_prediction = cross_val_predict(clf, X, y, cv=50, method='predict_proba', n_jobs=-1)
			cv_probabilities.append(clf_prediction)
			index = index + 1
		indx = 1
		for model in self.base_models:
			print("INSIDE MAIN TRAIN")
			print(f'---------------------------fitting base model {indx}-------------------------')
			model.fit(X, y)
			indx = indx + 1
		cv_probabilities = np.hstack(cv_probabilities)
		return cv_probabilities


	def transform_base_models(self, X):
		probabilities = []
		for model in self.base_models:
			pred_classifier = get_probabilities(model, X)
			probabilities.append(pred_classifier)
		return np.hstack(probabilities)


	def create_x_blending(self, X,  X_meta_probabilities):
		X_blending = np.hstack((X, X_meta_probabilities))
		return X_blending
		

	def fit(self, X, y,node,index):	
		X_meta_probabilities = self.fit_base_models(X, y)
		X = self.create_x_blending(X, X_meta_probabilities)
		self.meta_classifier.fit(X,y)
		self.classes_ = np.unique(y)	
		return self


	def create_x_test_blending(self, X_test):
		x_pred_test = self.transform_base_models(X_test)
		X_test_blending = np.hstack((X_test, x_pred_test))
		return X_test_blending

	def transform_meta_classifier_probabilities(self, X_blend):	
		model_probabilities = get_probabilities(self.meta_classifier, X_blend)
		return model_probabilities

	def transform_meta_classifier_predictions(self, X_blend):
		model_predictions = get_predictions(self.meta_classifier, X_blend)
		return model_predictions

	def predict_proba(self, X_test):
		X_blend = self.create_x_test_blending(X_test)
		probabilities = self.transform_meta_classifier_probabilities(X_blend)
		return probabilities

	def predict(self, X_test):
		X_blend = self.create_x_test_blending(X_test)
		predictions = self.transform_meta_classifier_predictions(X_blend)
		return predictions
