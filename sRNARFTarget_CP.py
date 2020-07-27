
	#!pip install git+https://github.com/ModelOriented/pyCeterisParibus

	# This approach calculates individual_variable_profile for a given observation. For Example - if grid_point parameter is 
	# given the value of 100, it changes the value of 1 variable 100 times, keeping all other features value constant and gets the model prediction.
	#In plot it shows those prediction value for number of grid_points times.
	# S0 if total number of variables are 29 and grid_points was passed with a value of 100, 
	#it will calculate 2900 profiles (100 for each variable). Plots will show how prediction changes as a value of particular variable changes.
import sys
if((len(sys.argv) - 1) == 3 ):
	import shap
	import pickle
	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.ensemble import RandomForestClassifier
	from ceteris_paribus.explainer import explain
	from ceteris_paribus.profiles import individual_variable_profile
	from ceteris_paribus.profiles import CeterisParibus
	from ceteris_paribus.plots.plots import plot, plot_notebook

	RFModel = pickle.load(open('./PickledModelData/RFModel/sRNARFTargetModel.pickle', 'rb'))

	featureData = pd.read_csv("./sRNARFTargetResult/FeatureFile.csv", sep='\t')
	datarow = featureData[(featureData['sRNA_ID'] == sys.argv[1]) & (featureData['mRNA_ID'] == sys.argv[2])]
	data_for_prediction = datarow.iloc[:, 2:]

	trainX = pickle.load(open('./PickledModelData/RFData/trainX_sRNARFTarget.pkl', 'rb'))
	trainY = pickle.load(open('./PickledModelData/RFData/trainY_sRNARFTarget.pkl', 'rb'))

	data = np.array(trainX)
	yt = np.array(trainY)
	labels = yt.ravel()
	variable_names = data_for_prediction.columns

	predict_function = lambda X: RFModel.predict_proba(X)[::, 1]
	explainer_rf = explain(RFModel, variable_names, data, y = labels, predict_function=predict_function, label = "sRNARFTarget")

	#cp_profile = individual_variable_profile(explainer_rf, data_for_prediction, y = 1, grid_points = 100)
	cp_profile = individual_variable_profile(explainer_rf, data_for_prediction, grid_points = 200, variables = [sys.argv[3]])
	plot(cp_profile, show_profiles=True, show_residuals = True, show_rugs=True, height=700, width=750, yaxis_title='Prediction probablity for class 1', plot_title = 'Ceteris paribus profiles of feature '+ sys.argv[3] +' for '+ sys.argv[1] +'-'+ sys.argv[2] + ' pair interaction',  
		color='blue', size=3, alpha=0.5,  color_residuals='red', size_residuals=20, alpha_residuals=20, print_observations = True)


elif((len(sys.argv) - 1) < 3):

	print("Error: Required parameters not passed! Please pass all three parameters, sRNA ID, mRNA ID, variable name.")


else:

	print("Error: Only three parameter can be passed, sRNA ID, mRNA ID, variable name.")