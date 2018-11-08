"""
dist_adj_weight_decision.py

Predicted number of synaptic connections for each cell compared 
to the actual number. Predictions made using a logistic regression 
classifier model. Red line indicates perfect agreement between predicted 
and actual values. The residual is the distance from the data point to 
the line. Colors indicate the probability of observing a residual
as large or larger. p adj is a representative probability for all data 
points, computed using multiple hypothesis testing.

created: Christopher Brittin
date: 01 November 2018

"""
import sys
sys.path.append(r'./volumetric_analysis')
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import digamma

import matplotlib as mpl
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import train_test_split
from sklearn import metrics
from sklearn.cross_validation import cross_val_score


from connectome.load import from_db
from networks.stats import get_corresponding_edge_attr
from models.mutual_info import *
from figures.stats import plot_adj_syn_mi
from models.mht import *

mpl.rcParams['xtick.labelsize'] = 32
mpl.rcParams['ytick.labelsize'] = 32

db = 'N2U'
remove = ['VC01','VD01','VB01','VB02']
SCALE = 5*90*(1e-6)
KMAX = 100
THETA = 0

def get_data(G1,G2):
    #Get edge weights
    N = G2.ecount()
    data = np.zeros((N,2))
    for i in range(N):
        e = G2.es[i]
        data[i,0] = e['weight']
        if G1.are_connected(e.source,e.target):
            w = G1.es[G1.get_eid(e.source,e.target)]['weight']
            if w > THETA:            
                data[i,1] = 1
    data[:,0] *= SCALE
    data[:,0] = np.log(data[:,0])
    
    return data

def result_data(G1,G2,model):
    N = G2.vcount()
    data = np.zeros((N,4))
    for i in range(N):
        size = G2.degree(i)
        #actual = C.C.degree(i)
        actual = 0
        for e in G1.incident(i):
            if G1.es[e]['weight'] > THETA:
                actual += 1
        p = actual / float(size)
        var = size*p*(1-p)
        w = np.log(np.array(G2.es[G2.incident(i)]['weight']).reshape(-1,1)*SCALE)
        predict = np.sum(model.predict(w))
        #w[w >= 1.28] = 1
        #w[w < 1] = 0
        #predict = np.sum(w)
        data[i,:] = [size,actual,predict,var]
            
    
    return data

def run(fout=None):
    C = from_db(db,adjacency=True,chemical=True,electrical=True,remove=remove)
    C.C.reduce_to(C.A)
    C.E.reduce_to(C.A)
    N = C.C.ecount()

    C.C.to_undirected(combine_edges=sum)

    data = get_data(C.C,C.A)
    data = data[data[:,0].argsort()]

    X = data[:,0].reshape(-1,1)
    y = np.ravel(data[:,1])
    _x = np.linspace(-4,4,81).reshape(-1,1)

    # instantiate a logistic regression model, and fit with X and y
    model = LogisticRegression()
    model = model.fit(X, y)

    # check the accuracy on the training set
    print(model.score(X, y))
    print(y.mean())
    #print(model.predict_proba(_x))

    _data = result_data(C.C,C.A,model)

    """
    plt.figure()
    plt.plot(data[:,1],data[:,2],'bo')
    plt.plot([0,50],[0,50],'r-',linewidth=3)
    plt.xlim([0,50])
    plt.ylim([0,50])
    plt.show()
    """

    fig,ax = plt.subplots(1,1,figsize=(15,10))
    plot_actual_vs_predict(ax,_data,colorbar=True)
    ax.set_xlim([0,50])
    ax.set_ylim([0,50])
    ax.set_title('Predicted number of synaptic connections per cell',
                 fontsize=32,y=1.04)
    ax.set_xlabel('# actual connections',fontsize=32)
    ax.set_ylabel('# predicted connections',fontsize=32)
    plt.tight_layout()
    if fout: plt.savefig(fout)
    plt.show()

    """
    # evaluate the model by splitting into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3,random_state=0)
    model2 = LogisticRegression()
    model2.fit(X_train, y_train)

    # predict class labels for the test set
    predicted = model2.predict(X_test)
    print(predicted)

    # generate class probabilities
    probs = model2.predict_proba(X_test)
    print(probs)
    
    # generate evaluation metrics
    print(metrics.accuracy_score(y_test, predicted))
    print(metrics.roc_auc_score(y_test, probs[:, 1]))

    print(metrics.confusion_matrix(y_test, predicted))
    print(metrics.classification_report(y_test, predicted))
    
    # evaluate the model using 10-fold cross-validation
    scores = cross_val_score(LogisticRegression(), X, y, scoring='accuracy', cv=10)
    print(scores)
    print(scores.mean())
    """

if __name__ == "__main__":
    run()
