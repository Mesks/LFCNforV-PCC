#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
from random import shuffle
import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Activation, Dropout, Dense, Flatten, Input
from tensorflow.keras.models import Model
from sklearn.utils import class_weight
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import concatenate
from keras.utils.vis_utils import plot_model
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.random
import argparse
import locale
import os

seed = 246

# model-compile parameter sets
model_metrics = 'acc'
epochs = 500
batchs = 64
splits = 0.2
lr        = 3e-4
input_dim = 3
opt = Adam(learning_rate=lr,decay=3e-4/200)

concatenated_df=pd.read_csv("oriI_extraFeatures_Geo.csv")
XY = concatenated_df.values
for i in range(10):
    np.random.shuffle(XY)
X = XY[:,[0,2,3]]
Y = XY[:,[5]]
x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=splits, random_state=seed)

model = Sequential()
inputShape=(input_dim,)
model.add(Input(shape=inputShape))
x = Dense(10,activation="relu", kernel_initializer="RandomNormal", bias_initializer="RandomNormal")(model.output)
x = Dense(5,activation ="sigmoid", kernel_initializer="RandomNormal", bias_initializer="RandomNormal")(x)
# x = Dense(5,activation ="sigmoid", kernel_initializer="RandomNormal", bias_initializer="RandomNormal")(x)
x = Dense(1,activation ="sigmoid", kernel_initializer="RandomNormal", bias_initializer="RandomNormal")(x)
model = Model(inputs=[model.input],outputs=x)
model.compile(loss="mse",optimizer=opt,metrics=['acc'])

y_train = y_train.flatten()
class_weights = class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)
class_weights = dict(zip(np.unique(y_train),class_weights))

plot_model(model,to_file='FeaturesPlots/model.png',show_shapes=True)


# In[2]:


history = model.fit(x=[x_train],y=y_train, validation_data=([x_test], y_test), 
                    epochs=epochs, batch_size=batchs, class_weight=class_weights)

model.save_weights(r'weightANDlearningcurve/GeoIModule_model_3layer.h5')
eval_model=[]
eval_model.append(model.evaluate([x_test], y_test)[1])
print("\nTest Accuracy: %.4f" % eval_model[0])


# In[3]:


plt.plot(history.history['loss'],color='r')
plt.plot(history.history['val_loss'],color='g')
plt.plot(history.history['acc'],color='b')
plt.plot(history.history['val_acc'],color='k')
plt.title('T1: I Frame module learning curve (Geometry)')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train_loss', 'test_loss','train_acc', 'test_acc'], loc='upper left',bbox_to_anchor=(0,-0.3))
plt.savefig('FeaturesPlots/I_GeoTrainingCurve.jpg', bbox_inches='tight', dpi=1280)
plt.show()

import pickle
with open('weightANDlearningcurve/GeoIModule_history.txt', 'wb') as file_txt:
    pickle.dump(history.history, file_txt)


# In[4]:


np.set_printoptions(suppress=True)

g_weight1=model.get_weights()[0]
g_bias1=model.get_weights()[1]
g_weight2=model.get_weights()[2]
g_bias2=model.get_weights()[3]
g_weight3=model.get_weights()[4]
g_bias3=model.get_weights()[5]


print("\ng_weight1: ")
for a in g_weight1:
    for b in a:
        print(b,end=",")
        
print("\n\ng_bias1: ")
for a in g_bias1:
        print(a,end=",")
        
print("\n\ng_weight2: ")
for a in g_weight2:
    for b in a:
        print(b,end=",")
        
print("\n\ng_bias2: ")
for a in g_bias2:
        print(a,end=",")

print("\n\na_weight3: ")
for a in g_weight3:
    for b in a:
        print(b,end=",")
        
print("\n\na_bias3: ")
for a in g_bias3:
        print(a,end=",")
        
# a_weight1=model.get_layer(index=0).get_weights()
# a_weight2=model.get_layer(index=1).get_weights()
        
# print(a_weight1)
# print(a_weight2)

