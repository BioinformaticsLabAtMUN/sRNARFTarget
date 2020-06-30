#!/usr/bin/env python

import pandas as pd
import numpy  as np 
from scipy import interp
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import precision_recall_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split
from IPython.display import display
from sklearn.model_selection import cross_validate
from sklearn.model_selection import ShuffleSplit
from sklearn.metrics import r2_score
from collections import defaultdict
from sklearn.model_selection import KFold

data1 = pd.read_csv('./Data/Model_Training/Tri_difference_distance_748class_shufffled_Nodistance.txt', sep='\t', header=None)

df1 = pd.DataFrame(data = data1)
X_features = df1.iloc[:, :-1].values #All columns except  last
Y_label = df1.iloc[:,  -1].values #Only last column

#Find best KNN parameters - model
def findbestKNNparameter():
    knn = KNeighborsClassifier()
    k_range = list(range(1, 50))
    param_grid = { 
    'n_neighbors': k_range,
    'weights': ['distance']
   }
    grid = GridSearchCV(knn, param_grid, cv=10, scoring='accuracy',n_jobs = -1,iid = False)
    grid.fit(X_features, Y_label)
    print(grid.best_score_)
    print(grid.best_params_)
    print(grid.best_estimator_)
    return grid.best_estimator_

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
avgpre = []
aucpr_vs_t=[] 
aucpr_vs_t1=[]
knnclf = findbestKNNparameter()


skfcv = StratifiedKFold(n_splits=10)
print(skfcv.get_n_splits(X_features, Y_label))
i = 0
for train_index, test_index in skfcv.split(X_features,Y_label):
    X_train, X_test, y_train, y_test = X_features[train_index], X_features[test_index],Y_label[train_index], Y_label[test_index]
    knnclf.fit(X_train, y_train)
    probas_ = knnclf.predict_proba(X_test)
    y_pred = knnclf.predict(X_test)
    
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
  
    precision, recall, _ = precision_recall_curve(y_test, probas_[:, 1])
    fpr, tpr, _          = roc_curve(y_test,  probas_[:, 1])
    roc_auc              = roc_auc_score(y_test,  probas_[:, 1])
    ave_precision        = average_precision_score(y_test,  probas_[:, 1])
    avgpre.append(ave_precision)
    i += 1
    
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',label='Random', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('k-nearest neighbors Triucleotides difference')

plt.legend(loc="upper left")
plt.show()

plt.savefig('KNN_ROC.png')

aucpr_vs_t = np.mean(avgpre)
aucpr_vs_t1 = np.std(avgpre)
print('mean Average Precision=%.2f' % (aucpr_vs_t))
print('std Average Precision=%.2f' % (aucpr_vs_t1))

