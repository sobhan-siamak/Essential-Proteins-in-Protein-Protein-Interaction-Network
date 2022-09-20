


# @copy by sobhan siamak

import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt
import time
from datetime import datetime
import networkx as nx
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report
# from sklearn import cross_validation
from sklearn.model_selection import cross_val_score,KFold, cross_validate
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
# import scikitplot as skplt

start_time = datetime.now()

db = pd.read_csv("ppi_dmelanogaster.csv")
dbmain = pd.read_csv("comppi--compartments--tax_dmelanogaster_loc_all.txt", sep='\t')
labels = pd.read_csv("Labels.csv")
labels = labels[['locus','essentiality status']]
print(labels.head())

# la1 = pd.DataFrame(la)
# lb1 = pd.DataFrame(lb)
# la1.to_csv("goa.csv" , index=False)
# lb1.to_csv("gob.csv", index=False)
goa = dbmain.iloc[0:,3]
gob = dbmain.iloc[0:,11]

la = [x.split("|") for x in goa]
lb = [x.split("|") for x in gob]

alpha = 1 #coif of gene ontology
beta = 1 #coif of nieghborhood
gama = 0.001 #coif of score

print("la is:", la[0])
# print(goa)
print("lb is :", lb[0])
intersect = np.intersect1d(la[0],lb[0])
print(intersect)

go = []
for i in range(len(la)):
  go.append(len(np.intersect1d(la[i], lb[i]))/ len(la[i]+lb[i]))


######Change weigth #1

print(db.head())
db['go'] = go
db['Interaction Score'] = gama*db['Interaction Score'] + alpha*db['go']
db = db.drop(['go'], axis=1)
print(db.head())

db2 = np.array(db)
Pra = db.iloc[1:,0]
Prb = db.iloc[1:,1]
Sc = db.iloc[1:,2]
prota = db2[0:,0]
protb = db2[0:,1]
score = db2[0:,2]
Pa = pd.read_csv("ProtA.tab", sep='\t')
Pb = pd.read_csv("ProtB.tab", sep='\t')
frames =[Pa,Pb]
proteins = (pd.concat(frames))
proteins = proteins.drop_duplicates()
print(proteins.head())




#Labeling Dataset
prlabels = []
prot = proteins.iloc[:,1]

for i in prot:
    if(sum(labels['locus'].isin([i])) == 0):
        prlabels.append(['unknown'])
        continue
    temp = labels[labels['locus'].isin([i])]
    lbl = list(temp['essentiality status'])
    prlabels.append(lbl)

for i in range(len(prlabels)):
    if prlabels[i] == ['NE,NE']:
        prlabels[i] = ['NE']
    if prlabels[i] == ['E,NE']:
        prlabels[i] = ['E']
    if prlabels[i] == ['NE,E']:
        prlabels[i] = ['E']
    if prlabels[i] == ['E,E']:
        prlabels[i] = ['E']
    if prlabels[i] == []:
        prlabels[i] = ['unknown']

print(prlabels)
print(len(prlabels))
proteins['Labels'] = np.squeeze(prlabels)
print(proteins.head())
FProteins = proteins.drop('To', axis=1)
FProteins = FProteins.rename(columns={"From": "Proteins"})
FProteins = FProteins.drop_duplicates()
print(FProteins.head())
print(FProteins)
fpcsv = pd.DataFrame(FProteins)
fpcsv.to_csv("FPR.csv", index=False)

# unknown = sum(proteins['Labels'].isin(['unknown']))
# essential = sum(proteins['Labels'].isin(['E']))
# nessential = sum(proteins['Labels'].isin(['NE']))

unknown = sum(FProteins['Labels'].isin(['unknown']))
essential = sum(FProteins['Labels'].isin(['E']))
nessential = sum(FProteins['Labels'].isin(['NE']))


print("the number of unknown of proteins are:", unknown)
print("the number of essential of proteins are:", essential)
print("the number of not essential of proteins are:", nessential)




#Creation Graph
G2 = nx.from_pandas_edgelist(db, "Interactor A", "Interactor B", ['Interaction Score'], create_using= nx.Graph())
print(G2.nodes())
print(len(G2.nodes()))


nodes = list(G2.nodes())
nodelabels = []
for i in nodes:
    tmp = FProteins[FProteins['Proteins'].isin([i])]
    lbls = list(tmp['Labels'])
    nodelabels.append(lbls)

print("length nodelabels is:", len(nodelabels))
print(G2.nodes())
print(len(G2.nodes()))
Fnodelabels = []
for i in range(len(nodelabels)):
    if (nodelabels[i] == ['unknown']):
        G2.remove_node(nodes[i])
    else:
        Fnodelabels.append(nodelabels[i])
print(G2.nodes())


FNlabels =[]
for i in Fnodelabels:
    if i == ['NE']:
        FNlabels.append(0)
    else:FNlabels.append(1)
print(FNlabels)
print(len(FNlabels))

print("length Fnodelabels is:", len(Fnodelabels))
print("length G2.nodes is:", len(G2.nodes()))
Fnodes = G2.nodes()




######Centrality
central1 = nx.eigenvector_centrality(G2,weight = 'Interaction Score', max_iter=1000)
central1 = list(central1.values())
# central11 = pd.DataFrame(central1)
# central11.to_csv("central1.csv", index=False)
print(central1)
central2 = nx.degree_centrality(G2)
central2 = list(central2.values())
# central22 = pd.DataFrame(central2)
# central22.to_csv("central2.csv", index=False)
print(central2)
# central2 = nx.degree_pearson_correlation_coefficient(G2,x='P26802', y='Q9W5A9')
# print("pearson corralation is:", central2)
# print(type(central1))
#####High Time Complexity
# central3 = nx.closeness_centrality(G2)
# central3 = list(central3.values())
# central33 = pd.DataFrame(central3)
# central33.to_csv("central3.csv", index=False)
# print(central3)

#####High Time Complexity
# central4 = nx.betweenness_centrality(G2)
# central4 = list(central4.values())
# central44 = pd.DataFrame(central4)
# central44.to_csv("central4.csv", index=False)
# print(central4)

#####High Time Complexity
# central5 = nx.subgraph_centrality(G2)
# central5 = list(central5.values())
# central55 = pd.DataFrame(central5)
# central55.to_csv("central5.csv", index=False)
# print(central5)


###########Classifiers
# DecisionTreeClassifier(max_depth=5),\
# RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),


# central7 = nx.clustering(G2)
# print("clustering is :", central7)



print(len(Fnodes))
print(len(central1))
print(len(Fnodelabels))
FinalP = pd.DataFrame()
FinalP['Proteins'] = Fnodes
FinalP['C1'] = central1
FinalP['C2'] = central2
c3 = pd.read_csv("central3.csv")
c4 = pd.read_csv("central4.csv")
c5 = pd.read_csv("central5.csv")
print(c3.head())
print(c4.head())
print(c5.head())
FinalP['C3'] = c3
FinalP['C4'] = c4
FinalP['C5'] = c5
#labeling

FinalP['Labels'] = FNlabels
print(FinalP.head())
FinalP.to_csv("Finalp.csv", index=False)



#########################################################################
#Classification

FPro = pd.read_csv("Finalp.csv", header=None)
print(FPro.head())
x = FPro.iloc[1:,1:-1]
y = FPro.iloc[1:,-1]
target_names = ['0', '1']
# print(x.iloc[:,0])
# print(y)
# print(len(y))
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.20)
# clf = RandomForestClassifier(n_estimators=1000, max_depth=4,random_state=5,class_weight ='imbalanced')
# clf.fit(X_train, y_train)
# result = clf.predict(X_test)
# print(classification_report(y_test, result,target_names=target_names))


clf = KNeighborsClassifier(n_neighbors=5)
clf.fit(X_train, y_train)
result2 = clf.predict(X_test)
print(classification_report(y_test, result2, target_names=target_names))

clf1 = SVC(gamma=2, C=1)
# clf1 = MLPClassifier(alpha=1, max_iter=1000)
# clf1 = GaussianNB()
clf1.fit(X_train, y_train)
result3 = clf1.predict(X_test)
print(classification_report(y_test, result3, target_names=target_names))

clf2 = LogisticRegression()
clf2.fit(X_train,y_train)
result4 = clf2.predict(X_test)
print(classification_report(y_test, result4, target_names=target_names))

clf4 = MLPClassifier(alpha=0.1, max_iter=1000)
clf4.fit(X_train,y_train)
result5 = clf4.predict(X_test)
print(classification_report(y_test, result5, target_names=target_names))

clf5 = GaussianNB()
clf5.fit(X_train,y_train)
result6 = clf5.predict(X_test)
print(classification_report(y_test, result6, target_names=target_names))



# score = cross_val_score(clf, x, y, cv=5)
# print(score)
######KNN
scoring = ['precision_macro', 'recall_macro', 'accuracy']
scores2 = cross_validate(clf, x, y, cv=5, scoring=scoring)
print("recall for 5-fold cross validation is:")
print(scores2['test_recall_macro'])
yknn1 = scores2['test_recall_macro']
print("precision for 5-fold cross validation is:")
print(scores2['test_precision_macro'])
yknn2 = scores2['test_precision_macro']
print("accuracy for 5-fold cross validation is:")
print(scores2['test_accuracy'])
yknn3 = scores2['test_accuracy']


print("++++++++++++++++++++++++++++++++++++++++++")
##### Gaussian Naive Bayse
scoring = ['precision_macro', 'recall_macro', 'accuracy']
scores3 = cross_validate(clf5, x, y, cv=5, scoring=scoring)
print("recall for 5-fold cross validation in GNB is:")
print(scores3['test_recall_macro'])
ygnb1 = scores3['test_recall_macro']
print("precision for 5-fold cross validation in GNB is:")
print(scores3['test_precision_macro'])
ygnb2 = scores3['test_precision_macro']
print("accuracy for 5-fold cross validationin GNB is:")
print(scores3['test_accuracy'])
ygnb3 = scores3['test_accuracy']




#Plot
#5-Fold Curve
xfold = [1,2,3,4,5]
plt.plot(xfold, yknn1, label='Recall')
plt.plot(xfold, yknn2, label='Precision')
plt.plot(xfold, yknn3, label='Accuracy')
plt.xlabel('5-Fold')
plt.ylabel('Measurment')
plt.title('KNN 5-Fold Curve')
plt.legend(loc='best')
plt.show()

xfold = [1,2,3,4,5]
plt.plot(xfold, ygnb1, label='Recall')
plt.plot(xfold, ygnb2, label='Precision')
plt.plot(xfold, ygnb3, label='Accuracy')
plt.xlabel('5-Fold')
plt.ylabel('Measurment')
plt.title(' Guassian Naive Bayse 5-Fold Curve')
plt.legend(loc='best')
plt.show()

#####ROC Curve
# fpr_knn_lm, tpr_knn_lm = roc_curve(y_test, result2)
# fpr_svc_lm, tpr_svc_lm = roc_curve(y_test, result3)
# fpr_lr_lm, tpr_lr_lm = roc_curve(y_test, result4)
#
# plt.figure(1)
# plt.plot([0, 1], [0, 1], 'k--')
# plt.plot(fpr_knn_lm, tpr_knn_lm, label='KNN')
# plt.plot(fpr_svc_lm, tpr_svc_lm, label='SVC')
# plt.xlabel('False positive rate')
# plt.ylabel('True positive rate')
# plt.title('ROC curve')
# plt.legend(loc='best')
# plt.show()





end_time = datetime.now()
print('Time for running this Code is: {}'.format(end_time - start_time))