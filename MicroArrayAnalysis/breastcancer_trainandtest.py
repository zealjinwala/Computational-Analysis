# Author: Zeal Jinwala
# Date: Sunday, May 08, 2022

# Purpose: evaluation function that trains an SVM 
# using Xtrain & Ttrain, where Xtrain contains gene expression data for a subset of the samples, 
# and Ttrain is a binary vector of class labels (indicating cancer relapse status) and 
# calculates the number of errors on the test data Xtest & Ttest. 

from sklearn import svm
def hwmaml_breastcancer_trainandtest(xTrain,yTrain,Xtest,Ttest):
    clf = svm.SVC(kernel='rbf') 
    clf.fit(xTrain, yTrain)
    y_pred = clf.predict(Xtest)
    numerrors = sum(y_pred != Ttest)
    return numerrors