#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:42:42 2019

@author: ekrem
"""

from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from keras.layers import GlobalMaxPooling2D
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import BatchNormalization
from keras.layers import Activation
from keras.layers import GaussianNoise
from keras.optimizers import Adam, SGD
from keras.layers import concatenate
from keras import Input, Model
from keras.regularizers import l2

'''
When only the consensus repeat is available
'''
def build_parallel_classifier_R(seq_length, kernel_width = [4, 6, 8, 12, 16]):
    
    seq_height = 4
    reg = 0.05
    #noise_stddev = 0.015
    num_feature_maps = 32
    
    main_input = Input(shape=(seq_height, seq_length, 1))
    
    ''''''
    layer_1p = Conv2D(num_feature_maps, (seq_height, kernel_width[0]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_1p = GaussianNoise(noise_stddev)(layer_1p)
    layer_1p = Activation('relu')(layer_1p)
    layer_1p = BatchNormalization()(layer_1p)
    pool_1p = GlobalMaxPooling2D()(layer_1p)
    
    layer_2p = Conv2D(num_feature_maps, (seq_height, kernel_width[1]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_2p = GaussianNoise(noise_stddev)(layer_2p)
    layer_2p = Activation('relu')(layer_2p)
    layer_2p = BatchNormalization()(layer_2p)
    pool_2p = GlobalMaxPooling2D()(layer_2p)
    
    layer_3p = Conv2D(num_feature_maps, (seq_height, kernel_width[2]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_3p = GaussianNoise(noise_stddev)(layer_3p)
    layer_3p = Activation('relu')(layer_3p)
    layer_3p = BatchNormalization()(layer_3p)
    pool_3p = GlobalMaxPooling2D()(layer_3p)
    
    layer_4p = Conv2D(num_feature_maps, (seq_height, kernel_width[3]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_4p = GaussianNoise(noise_stddev)(layer_4p)
    layer_4p = Activation('relu')(layer_4p)
    layer_4p = BatchNormalization()(layer_4p)
    pool_4p = GlobalMaxPooling2D()(layer_4p)
    
    layer_5p = Conv2D(num_feature_maps, (seq_height, kernel_width[4]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_5p = GaussianNoise(noise_stddev)(layer_5p)
    layer_5p = Activation('relu')(layer_5p)
    layer_5p = BatchNormalization()(layer_5p)
    pool_5p = GlobalMaxPooling2D()(layer_5p)
    
    concatenated = concatenate([pool_1p, pool_2p, pool_3p, pool_4p, pool_5p], axis = 1, name = 'cutoff_layer')
    
    x = Dense(256, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(concatenated)
    x = Activation('relu')(x)
    x = Dropout(rate = 0.5)(x)
    
    x = Dense(32, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(x)
    x = Activation('relu')(x)
    x = Dropout(rate = 0.5)(x)
    
    xs = Dense(units = 1, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(x)
    out = Activation('sigmoid')(xs)

    classifier = Model(main_input, out)
    
    optim = SGD(decay=1e-4)
    classifier.compile(optimizer = optim, loss = 'binary_crossentropy', metrics = ['accuracy'])
    return classifier

'''
When a sequence is available
'''
def build_parallel_classifier_A(seq_length, kernel_width = [4, 6, 8, 12, 16]):
    
    seq_height = 5
    reg = 0.05
    #noise_stddev = 0.015
    num_feature_maps = 32
    
    main_input = Input(shape=(seq_height, seq_length, 1))
    side_input = Input(shape=(7,))
    
    ''''''
    layer_1p = Conv2D(num_feature_maps, (seq_height, kernel_width[0]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_1p = GaussianNoise(noise_stddev)(layer_1p)
    layer_1p = Activation('relu')(layer_1p)
    layer_1p = BatchNormalization()(layer_1p)
    pool_1p = GlobalMaxPooling2D()(layer_1p)
    
    layer_2p = Conv2D(num_feature_maps, (seq_height, kernel_width[1]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_2p = GaussianNoise(noise_stddev)(layer_2p)
    layer_2p = Activation('relu')(layer_2p)
    layer_2p = BatchNormalization()(layer_2p)
    pool_2p = GlobalMaxPooling2D()(layer_2p)
    
    layer_3p = Conv2D(num_feature_maps, (seq_height, kernel_width[2]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_3p = GaussianNoise(noise_stddev)(layer_3p)
    layer_3p = Activation('relu')(layer_3p)
    layer_3p = BatchNormalization()(layer_3p)
    pool_3p = GlobalMaxPooling2D()(layer_3p)
    
    layer_4p = Conv2D(num_feature_maps, (seq_height, kernel_width[3]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_4p = GaussianNoise(noise_stddev)(layer_4p)
    layer_4p = Activation('relu')(layer_4p)
    layer_4p = BatchNormalization()(layer_4p)
    pool_4p = GlobalMaxPooling2D()(layer_4p)
    
    layer_5p = Conv2D(num_feature_maps, (seq_height, kernel_width[4]), padding = 'valid', use_bias = False, kernel_regularizer=l2(reg))(main_input)
    #layer_5p = GaussianNoise(noise_stddev)(layer_5p)
    layer_5p = Activation('relu')(layer_5p)
    layer_5p = BatchNormalization()(layer_5p)
    pool_5p = GlobalMaxPooling2D()(layer_5p)
    
    concatenated = concatenate([pool_1p, pool_2p, pool_3p, pool_4p, pool_5p], axis = 1, name = 'cutoff_layer')
    
    x = Dense(256, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(concatenated)
    x = Activation('relu')(x)
    x = Dropout(rate = 0.5)(x)
    
    x = Dense(32, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(x)
    x = Activation('relu')(x)
    x = Dropout(rate = 0.5)(x)
    
    s = Dense(16, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(side_input)
    s = Activation('relu')(s)
    s = Dropout(rate = 0.5)(s)
    
    xs = concatenate([x, s], axis = 1)
    xs = Dense(units = 1, kernel_regularizer=l2(reg), bias_regularizer=l2(reg))(xs)
    out = Activation('sigmoid')(xs)
    
    classifier = Model([main_input, side_input], out)
    
    optim = SGD(decay=1e-4)
    classifier.compile(optimizer = optim, loss = 'binary_crossentropy', metrics = ['accuracy'])
    return classifier