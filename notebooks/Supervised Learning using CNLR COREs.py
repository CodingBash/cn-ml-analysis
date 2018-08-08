
# coding: utf-8

# # Supervised Learning of Drug Response using CORES from Copy Number Log Ratio

# ### Import Python source code

# In[1]:


"""
Created on Thu Jul 26 12:21:38 2018

@author: bbece
"""

from __future__ import division, print_function, unicode_literals
import numpy as np
import os
from IPython.display import display, HTML

from pprint import pprint
np.random.seed(42)

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
            
import pandas as pd
import scipy

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import MinMaxScaler

from sklearn import decomposition

from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score

from sklearn.tree import export_graphviz

import math



# ### Load training set matrix

# In[43]:

labeled_matrix_training_set = pd.read_csv("../mlOutput/geneTrainingSet_8_5_2018_1.csv") # Cancer gene feature set
#labeled_matrix_training_set = pd.read_csv("../mlOutput/slicingTrainingSet_8_3_2018_1.csv") # Slicing features on A and D combined
#labeled_matrix_training_set = pd.read_csv("../mlOutput/coreTrainingSet_8_2_2018_2.csv") # Segment features on A and D combined
#labeled_matrix_training_set = pd.read_csv("../mlOutput/coreTrainingSet_8_2_2018_1.csv") # Segment features on A and D merged
#labeled_matrix_training_set.columns.values[0] = "sampleId"
labeled_matrix_training_set = labeled_matrix_training_set.drop([labeled_matrix_training_set.columns[0]], axis = 1)
labels = list(range(0,5))


# In[44]:

display(labeled_matrix_training_set.head(25))


# In[45]:

X = labeled_matrix_training_set.copy().drop(labeled_matrix_training_set.columns[labels], axis = 1)
y = labeled_matrix_training_set.copy()[labeled_matrix_training_set.columns[labels]]


# In[46]:

display(X.head())


# In[47]:

display(y.head(15))


# In[48]:

from sklearn.model_selection import train_test_split

all_X_TRAIN, all_X_TEST, all_Y_TRAIN, all_Y_TEST = train_test_split(X, y, test_size=0.20, random_state=42)
# TODO: train_test must be split on amount of NAs as well!


# In[49]:

display(all_X_TRAIN.head())
display(all_Y_TRAIN.head())


# In[50]:

display(all_X_TEST.head())
display(all_Y_TEST.head())


# ## Visualize ML Results

# In[51]:

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


# In[52]:

def retrieve_pipelines(model_name, ml_model, scaling_method = "standard"):
    
    labelScaler = StandardScaler() if scaling_method == "standard" else MinMaxScaler()
    featureScaler = StandardScaler() if scaling_method == "standard" else MinMaxScaler()
    
    if scaling_method == "none":
        labelPipeline = Pipeline([
         ('imputer', Imputer(axis=0,strategy="median")),
        ])


        featureSetPipeline = Pipeline([
            ('imputer', Imputer(axis=0,strategy="median")),
        ])
    else:
        labelPipeline = Pipeline([
         ('imputer', Imputer(axis=0,strategy="median")),
         ('scaler', labelScaler),
        ])


        featureSetPipeline = Pipeline([
            ('imputer', Imputer(axis=0,strategy="median")),
            ('scaler', featureScaler),
        ])
    
    modelPipeline = Pipeline([
        #('imputer', Imputer(axis=0,strategy="median")),
        #('scaler', StandardScaler()),
        #('scaler', MinMaxScaler()),
        #("pca", decomposition.PCA(n_components=10)),
        (model_name,  ml_model)
    ])
    
    return (modelPipeline, labelPipeline, featureSetPipeline)

def imputer_inverse_transform(pre_data, post_data):
    na_indices = np.where(np.isnan(pre_data))[0]
    post_data[na_indices] = float('NaN')
    return post_data
    
def remove_NAs(X, y, label):
    label_y = y[[y.columns[label]]]
    na_indices = label_y[label_y.columns[0]].index[label_y[label_y.columns[0]].apply(np.isnan)]
    y_nonNA = label_y.copy().drop(na_indices)
    X_nonNA = X.copy().drop(na_indices)
    return X_nonNA, y_nonNA
    
def train_and_test(labelPipeline, modelPipeline, X_TRAIN, X_TEST, this_y_train, this_y_test, scaling = True):
    modelPipeline.fit(X_TRAIN,this_y_train)

    y_prediction = modelPipeline.predict(X_TEST)
    
    if scaling ==  True:
        y_prediction = labelPipeline.named_steps['scaler'].inverse_transform(y_prediction)
    y_prediction = imputer_inverse_transform(this_y_test, y_prediction)

    if scaling == True:
        y_test_np = labelPipeline.named_steps['scaler'].inverse_transform(this_y_test)
    else:
        y_test_np = this_y_test
    y_test_np = y_test_np[~np.isnan(y_test_np)]
    y_prediction = y_prediction[~np.isnan(y_prediction)]
    return (y_test_np, y_prediction)


def simple_score(y_test_np, y_prediction):
    rmse = np.sqrt(mean_squared_error(y_test_np, y_prediction))
    r = scipy.stats.pearsonr(y_test_np, y_prediction)
    t = scipy.stats.spearmanr(y_test_np, y_prediction)
    return (rmse, r, t)
    
def visualize(y_test_np, y_prediction):
    plt.plot(y_test_np, y_prediction, 'bo')
    abline(1,0)
    plt.ylabel("Prediction")
    plt.xlabel("Label")
    plt.show()
    
def display_set(X_TRAIN, X_TEST, Y_TRAIN, Y_TEST):
    TRAIN = pd.concat([Y_TRAIN, X_TRAIN], axis=1)
    TEST = pd.concat([Y_TEST, X_TEST], axis=1)
    pprint("TRAIN")
    display(TRAIN)
    pprint("TEST")
    display(TEST)

def cv_score(XYpipeline, X_TRAIN, this_y_train_tr):
    scores = cross_val_score(XYpipeline, X_TRAIN, this_y_train_tr,
                             scoring = "neg_mean_squared_error", cv=10)
    scores = Ypipeline.named_steps['normalizer'].inverse_transform(scores)
    return scores


# ### Visualize ML results using Ridge Regression

# In[53]:

for label in labels:
    modelPipeline, labelPipeline, featureSetPipeline  = retrieve_pipelines("ridge_model", Ridge(alpha = 0.5), "normal")
    
    # Impute samples where label is NA
    X_nonNA, y_nonNA = remove_NAs(X, y, label)
    
    X_nonNA_tr = featureSetPipeline.fit_transform(X_nonNA)
    y_nonNA_tr = labelPipeline.fit_transform(y_nonNA)
    
    X_nonNA_tr = pd.DataFrame(data=X_nonNA_tr, index=X_nonNA.index.values, columns = X_nonNA.columns.values)
    y_nonNA_tr = pd.DataFrame(data=y_nonNA_tr, index=y_nonNA.index.values, columns = y_nonNA.columns.values)
    
    X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(X_nonNA_tr, y_nonNA_tr, test_size=0.10)
    
    
    display_set(X_TRAIN, X_TEST, Y_TRAIN, Y_TEST)

    y_test_np, y_prediction = train_and_test(labelPipeline, modelPipeline, X_TRAIN, X_TEST, Y_TRAIN.values, Y_TEST.values)

    rmse, r, t = simple_score(y_test_np, y_prediction)
    
    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))

    #scores = cv_score(XYpipeline, X_TRAIN, this_y_train_tr)
    
    #print("CV Scores: " + str(scores))
    #print("CV Mean: " + str(scores.mean()))
    #print("CV STD: " + str(scores.std()))
    print(y_prediction)
    visualize(y_test_np, y_prediction)


# ### Visualize ML results using Random Forest Regressor

# In[204]:

for label in labels:
    modelPipeline, labelPipeline, featureSetPipeline  = retrieve_pipelines("rfs_model", RandomForestRegressor(n_estimators=500, max_leaf_nodes=16, n_jobs=8), "normal")
    
    # Impute samples where label is NA
    X_nonNA, y_nonNA = remove_NAs(X, y, label)
    
    X_nonNA_tr = featureSetPipeline.fit_transform(X_nonNA)
    y_nonNA_tr = labelPipeline.fit_transform(y_nonNA)
    
    X_nonNA_tr = pd.DataFrame(data=X_nonNA_tr, index=X_nonNA.index.values, columns = X_nonNA.columns.values)
    y_nonNA_tr = pd.DataFrame(data=y_nonNA_tr, index=y_nonNA.index.values, columns = y_nonNA.columns.values)
    
    X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(X_nonNA_tr, y_nonNA_tr, test_size=0.15)
    
    
    display_set(X_TRAIN, X_TEST, Y_TRAIN, Y_TEST)

    y_test_np, y_prediction, this_y_train_tr = train_and_test(labelPipeline, XYpipeline, X_TRAIN, X_TEST, Y_TRAIN.values, Y_TEST.values)

    rmse, r, t = simple_score(y_test_np, y_prediction)
    
    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))

    #scores = cv_score(XYpipeline, X_TRAIN, this_y_train_tr)
    
    #print("CV Scores: " + str(scores))
    #print("CV Mean: " + str(scores.mean()))
    #print("CV STD: " + str(scores.std()))
    print(y_prediction)
    visualize(y_test_np, y_prediction)


# ### Bootstrap Regression Model

# In[57]:

for label in labels:
    modelPipeline, labelPipeline, featureSetPipeline  = retrieve_pipelines("ridge_model", Ridge(alpha = 0.5), "standard")
    
    # Impute samples where label is NA
    X_nonNA, y_nonNA = remove_NAs(X, y, label)
    
    X_nonNA_tr = featureSetPipeline.fit_transform(X_nonNA)
    y_nonNA_tr = labelPipeline.fit_transform(y_nonNA)
    
    X_nonNA_tr = pd.DataFrame(data=X_nonNA_tr, index=X_nonNA.index.values, columns = X_nonNA.columns.values)
    y_nonNA_tr = pd.DataFrame(data=y_nonNA_tr, index=y_nonNA.index.values, columns = y_nonNA.columns.values)

    
    all_y_test_np = np.array([])
    all_y_prediction = np.array([])
    for i in range(1,10):
        X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(X_nonNA_tr, y_nonNA_tr, test_size=0.05)

        y_test_np, y_prediction = train_and_test(labelPipeline, modelPipeline, X_TRAIN, X_TEST, Y_TRAIN.values, Y_TEST.values, scaling = False)

        all_y_test_np = np.append(all_y_test_np, y_test_np)
        all_y_prediction = np.append(all_y_prediction, y_prediction)
        
    rmse, r, t = simple_score(all_y_test_np, all_y_prediction)
    
    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))
    visualize(all_y_test_np, all_y_prediction)


# ### Bootstrap Random Forest Model

# In[40]:

for label in labels:
    X_nonNA, y_nonNA = remove_NAs(X, y, label)
    
    all_y_test_np = np.array([])
    all_y_prediction = np.array([])
    for i in range(1,10):
        X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(X_nonNA, y_nonNA, test_size=0.10)

        Ypipeline, XYpipeline = retrieve_pipelines("rfs_model", RandomForestRegressor(n_estimators=1000, max_leaf_nodes=8, n_jobs=4))

        y_test_np, y_prediction, _ = train_and_test(Ypipeline, XYpipeline, X_TRAIN, X_TEST, Y_TRAIN.values, Y_TEST.values)
        all_y_test_np = np.append(all_y_test_np, y_test_np)
        all_y_prediction = np.append(all_y_prediction, y_prediction)
        
    rmse, r, t = simple_score(all_y_test_np, all_y_prediction)
    
    print("RMSE: " + str(rmse))
    print("Pearson: " + str(r))
    print("Spearman: " + str(t))
    visualize(all_y_test_np, all_y_prediction)


# In[ ]:



