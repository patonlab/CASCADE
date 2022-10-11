""" Keras loss functions will choke on NaN inputs, which we'll often have for
unspecified inputs. These are just a couple simple loss functions that mask NaN
values in both the test and predicted tensors when computing the cost. """

import keras.backend as K
import tensorflow as tf

def masked_mean_squared_error(y_true, y_pred):
    mask = tf.is_finite(y_true)
    y_true_mask = tf.boolean_mask(y_true, mask)
    y_pred_mask = tf.boolean_mask(y_pred, mask)
    return K.mean(K.square(y_pred_mask - y_true_mask), axis=-1)


def masked_mean_absolute_error(y_true, y_pred):
    mask = tf.is_finite(y_true)
    y_true_mask = tf.boolean_mask(y_true, mask)
    y_pred_mask = tf.boolean_mask(y_pred, mask)
    return K.mean(K.abs(y_pred_mask - y_true_mask), axis=-1)
