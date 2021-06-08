#!/usr/bin/env python

#!pip install git+https://github.com/slundberg/shap.git

import sys
if((len(sys.argv) - 1) == 2 ):
        import shap
        import pickle
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from sklearn.ensemble import RandomForestClassifier
        
        # Load model and read data
        RFModel = pickle.load(open('./PickledModelData/RFModel/sRNARFTargetModel.pickle', 'rb'))
        
        data = pd.read_csv("./sRNARFTargetResult/FeatureFile.csv", sep='\t')
        
        datarow = data[(data['sRNA_ID'] == sys.argv[1]) & (data['mRNA_ID'] == sys.argv[2])]

        if (datarow.shape[0] == 0):
            print("No predictions found with", sys.argv[1], " and ", sys.argv[2])
            quit()
        
        data_for_prediction = datarow.iloc[:, 2:]
        
        #data_for_prediction_array = data_for_prediction.values.reshape(1, -1)
        #predicted_proba = RFModel.predict_proba(data_for_prediction)
        #RFprediction = RFModel.predict(data_for_prediction_array)
        #predict_proba = RFModel.predict_proba(data_for_prediction_array)
        
        import warnings
        warnings.filterwarnings("ignore")
        # Create SHAP explainer
        #explainer = shap.TreeExplainer(RFModel, data_for_prediction_array)
        explainer = shap.TreeExplainer(RFModel)
        
        # Get shap values for observation of interest
        shap_values = explainer.shap_values(data_for_prediction)
        
        #Decision plot
        expected_value = explainer.expected_value
        decisionhtml = shap.decision_plot(expected_value[1], shap_values[1], data_for_prediction.columns.tolist(), link='logit', show = False)
        plt.savefig('decisionPlot.pdf')
        plt.close()
        
        # Force plot
        forcehtml = shap.force_plot(explainer.expected_value[1], shap_values[1], data_for_prediction)
        shap.save_html(out_file = 'forcePlot.html', full_html=False, plot = forcehtml)
    
    
elif((len(sys.argv) - 1) < 2):
    print("Error: Required parameters not passed! Please pass two parameters, sRNA ID and mRNA ID.")
        
        
else:
    print("Error: Only two parameters can be passed. sRNA ID and mRNA ID.")
