# First keras analyses with in-silico flux analyses from previous project
import numpy as np
import pandas as pd
from tensorflow.contrib import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.linear_model import LassoLarsCV
from sklearn.preprocessing import LabelEncoder

# fix random seed for reproducibility
seed = 7
np.random.seed(7)

# Flux dataset from predicting bacterial growth conditions study
# import pdb; pdb.set_trace()
dataset1 = pd.read_csv("Data/syntheticFluxData.csv", delimiter=',', header=None)

# Later on build a function for processing dataset1
dataset = dataset1
all_input = dataset.iloc[:, 0:2381]
all_output = dataset.iloc[:, 2384]

# encode class values as integers (this needs to be convereted into funciton as well)
encoder = LabelEncoder()
encoder.fit(all_output)
encoded_Y = encoder.transform(all_output)
dummy_y = np_utils.to_categorical(encoded_Y)

# function to input number of input features (2382), dense layers(12) and outputs(1)
model = Sequential()
model.add(Dense(12, input_dim=2382, activation='relu'))
model.add(Dense(1, activation='softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

# batch size is how many samples you pass iteratively (efficient mem usage for low process computers)
estimator = KerasClassifier(model, epochs=200, batch_size=49)

# fold for cross-validation
kfold = KFold(n_splits=10, shuffle=True, random_state=seed)

# run model and use cv for metrics
X = np.array(all_input)
y = encoded_Y
# import pdb; pdb.set_trace()
results = cross_val_score(LassoLarsCV(), X, y, cv=kfold)
print(results)

#results = cross_val_score(estimator, all_input, dummy_y, cv=kfold)
