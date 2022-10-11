from keras.engine import Layer

from keras.layers import activations
from keras.layers import initializers
from keras.layers import regularizers
from keras.layers import constraints
from keras.layers import Lambda

import tensorflow as tf
import keras.backend as K

class MessageLayer(Layer):
    """ Implements the matrix multiplication message functions from Gilmer
    2017. This could probably be implemented as a series of other layers, but
    this is more convenient.

    """

    def __init__(self, dropout=0., reducer=None, **kwargs):
        """ 

        dropout : float between 0 and 1
            Whether to apply dropout to individual messages before they are
            reduced to each incoming atom.

        reducer : ['sum', 'mean', 'max', or 'min']
            How to collect incoming messages for each atom. In this library,
            I'm careful to only have messages be a function of the sending
            atom, so we can sort the connectivity matrix by recieving atom.
            That lets us use the `segment_*` methods from tensorflow, instead
            of the `unsorted_segment_*` methods.

        """

        self.dropout = dropout
        self.reducer = reducer

        reducer_dict = {
            None: tf.segment_sum,
            'sum': tf.segment_sum,
            'mean': tf.segment_mean,
            'max': tf.segment_max,
            'min': tf.segment_min
        }

        self._reducer = reducer_dict[reducer]

        super(MessageLayer, self).__init__(**kwargs)
    

    def call(self, inputs, training=None):
        """ Perform a single message passing step, returing the summed messages
        for each recieving atom.

        Inputs are [atom_matrix, bond_matrix, connectivity_matrix]

        atom_matrix : (num_atoms_in_batch, d)
            The input matrix of current hidden states for each atom

        bond_matrix : (num_bonds_in_batch, d, d)
            A matrix of current edge features, with each edge represented as a
            (dxd) matrix.

        connectivity : (num_bonds_in_batch, 2)
            A matrix of (a_i, a_j) pairs that indicates the bond in bond_matrix
            connecting atom_matrix[a_j] to atom_matrix[a_i].
            The first entry indicates the recieving atom.

        """

        atom_matrix, bond_matrix, connectivity = inputs

        # Gather the atom matrix so that each reciever node corresponds with
        # the given bond_matrix entry
        atom_gathered = tf.gather(atom_matrix, connectivity[:, 1])

        # Multiply the bond matrices by the gathered atom matrices
        messages = K.batch_dot(bond_matrix, atom_gathered)

        # Add dropout on a message-by-message basis if desired
        def add_dropout():
            if 0. < self.dropout < 1.:
                return K.dropout(messages, self.dropout)
            else:
                return messages

        dropout_messages = K.in_train_phase(
            add_dropout(), messages, training=training)

        # Sum each message along the (sorted) reciever nodes
        summed_message = self._reducer(dropout_messages, connectivity[:, 0])

        return summed_message


    def compute_output_shape(self, input_shape):
        """ Computes the shape of the output, which should be the same
        dimension as the first input, that atom hidden state """
        
        assert input_shape and len(input_shape) == 3
        assert input_shape[0][-1]  # atom hidden state dimension must be specified
        return input_shape[0]


    def get_config(self):
        config = {
            'dropout': self.dropout,
            'reducer': self.reducer,
        }
        base_config = super(MessageLayer, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))


class GatherAtomToBond(Layer):
    """ Reshapes the atom matrix (num_atoms_in_batch, d) to the bond matrix
    (num_bonds_in_batch, d) by reindexing according to which atom is involved
    in each bond.

    index : 0 or 1
        whether to gather the sending atoms (1) or recieving atoms (0) for each
        bond.

    """

    def __init__(self, index, **kwargs):
        self.index = index
        super(GatherAtomToBond, self).__init__(**kwargs)
    
    def call(self, inputs):
        atom_matrix, connectivity = inputs
        return  tf.gather(atom_matrix, connectivity[:, self.index])

    def compute_output_shape(self, input_shape):
        """ Computes the shape of the output,
        which should be the shape of the atom matrix with the length
        of the bond matrix """
        
        assert input_shape and len(input_shape) == 2
        return input_shape[0]
        
    def get_config(self):
        config = {
            'index': self.index,
        }
        base_config = super(GatherAtomToBond, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))

class Reducer(Layer):
    """ Superclass for reducing methods. 
    
    reducer : ['sum', 'mean', 'max', or 'min']
        How to collect elements for each atom or molecule. In this library,
        I'm careful to only have messages be a function of the sending
        atom, so we can sort the connectivity matrix by recieving atom.
        That lets us use the `segment_*` methods from tensorflow, instead
        of the `unsorted_segment_*` methods.

    """

    def __init__(self, reducer=None, **kwargs):

        self.reducer = reducer

        reducer_dict = {
            None: tf.segment_sum,
            'sum': tf.segment_sum,
            'unsorted_sum': tf.unsorted_segment_sum,
            'mean': tf.segment_mean,
            'unsorted_mean': tf.unsorted_segment_mean,
            'max': tf.segment_max,
            'min': tf.segment_min
        }

        self._reducer = reducer_dict[reducer]

        super(Reducer, self).__init__(**kwargs)

    def get_config(self):
        config = {'reducer': self.reducer}
        base_config = super(Reducer, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))

    def compute_output_shape(self, input_shape):
        assert input_shape and len(input_shape) == 2
        # Output shape is (n_graphs, atom_dim)
        return input_shape[0]


class ReduceAtomToMol(Reducer):
    """ Sum over all atoms in each molecule.

    Inputs

    atom_matrix : (num_atoms_in_batch, d)
        atom hidden states for each atom in the batch

    node_graph_indices : (num_atoms_in_batch,)
        A scalar for each atom representing which molecule in the batch the
        atom belongs to. This is generated by the preprocessor class, and
        essentially looks like [0, 0, 0, 1, 1] for a batch with a 3 atom
        molecule and a 2 atom molecule.
    """

    def call(self, inputs):

        atom_matrix, node_graph_indices = inputs
        return self._reducer(atom_matrix, node_graph_indices)

class ReduceBondToPro(Reducer):
    """
    Sums over bonds acoording to bond_index to get target bond properties
    
    Inputs

    bond_matrix : (num_bonds_in_batch, d)
        bond hidden states for each bond in the batch

    bond_index : (bond_atoms_in_batch, )
        A scalar for each bond representing which number in the target property 
        the bond links to. This is generated by the preprocessor class, and
        essentially looks like [-1,-1,0,0,-1,-1,1....]
    """

    def call(self, inputs):
        
        bond_matrix, bond_index, n_pro = inputs
        num_segments = tf.reduce_sum(n_pro)
        return self._reducer(bond_matrix, bond_index, num_segments)

    def compute_output_shape(self, input_shape):
        assert input_shape and len(input_shape) == 3
        return input_shape[0]

class ReduceAtomToPro(Reducer):
    """
    Sums over atoms acoording to atom_index to get target atom properties
    
    Inputs

    atom_matrix : (num_atoms_in_batch, d)
        atom hidden states for each atom in the batch

    atom_index : (atom_atoms_in_batch, )
        A scalar for each atom representing which number in the target property 
        the atom links to. This is generated by the preprocessor class, and
        essentially looks like [-1,-1,0,0,-1,-1,1....]
    """

    def call(self, inputs):
        
        atom_matrix, atom_index, n_pro = inputs
        num_segments = tf.reduce_sum(n_pro)
        return self._reducer(atom_matrix, atom_index, num_segments)

    def compute_output_shape(self, input_shape):
        assert input_shape and len(input_shape) == 3
        return input_shape[0]

class ReduceBondToAtom(Reducer):

    """ Sums over the incoming messages from all sender atoms.

    Inputs: 
    
    bond_matrix : (num_bonds_in_batch, d)
        A matrix of messages coming from each sender atom; one row for each
        bond/edge.

    connectivity : (num_bonds_in_batch, 2)
        A matrix of (a_i, a_j) pairs that indicates the bond in bond_matrix
        connecting atom_matrix[a_j] to atom_matrix[a_i].
        The first entry indicates the recieving atom.

    Again, I'm careful to only have the messages be a function of the sending
    node, such that we can use sorted methods in performing the reduction.

    """

    def call(self, inputs):

        bond_matrix, connectivity = inputs
        return self._reducer(bond_matrix, connectivity[:, 0])


class Squeeze(Layer):
    """ Keras forces inputs to be a vector per entry, so this layer squeezes
    them to a single dimension.

    I.e., node_graph_indices will have shape (num_atoms_in_batch, 1), while its
    easier to work with a vector of shape (num_atoms_in_batch,)
    """

    def call(self, inputs):
        return K.squeeze(inputs, 1)

    def compute_output_shape(self, input_shape):
        return input_shape[:-1]


class Embedding2D(Layer):
    """ Keras typically wants to embed items as a single vector, while for the
    matrix multiplication method of Gilmer 2017 we need a matrix for each bond
    type. This just implements that fairly simple extension of the traditional
    embedding layer.
    """

    def __init__(self, input_dim, output_dim,
                 embeddings_initializer='uniform',
                 embeddings_regularizer=None,
                 embeddings_constraint=None,
                 **kwargs):
        
        super(Embedding2D, self).__init__(**kwargs)

        self.input_dim = input_dim
        self.output_dim = output_dim
        
        self.embeddings_initializer = initializers.get(embeddings_initializer)
        self.embeddings_regularizer = regularizers.get(embeddings_regularizer)
        self.embeddings_constraint = constraints.get(embeddings_constraint)
        
    def build(self, input_shape):

        self.embeddings = self.add_weight(
            shape=(self.input_dim, self.output_dim, self.output_dim),
            initializer=self.embeddings_initializer,
            name='bond_embedding_weights', 
            regularizer=self.embeddings_regularizer,
            constraint=self.embeddings_constraint)
        
        self.built = True
        
    def call(self, inputs):
        return tf.nn.embedding_lookup(self.embeddings, inputs)
        
    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim, self.output_dim)

    def get_config(self):
        config = {
            'input_dim': self.input_dim,
            'output_dim': self.output_dim,
            'embeddings_initializer':
            initializers.serialize(self.embeddings_initializer),
            'embeddings_regularizer':
            regularizers.serialize(self.embeddings_regularizer),
            'embeddings_constraint':
            constraints.serialize(self.embeddings_constraint),
        }
        base_config = super(Embedding2D, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))


class EdgeNetwork(Layer):
    """ A layer to embed (bond_type, distance) pairs as a NxN matrix. 
    
    Inputs:
    units : dimension of the output matrix
    bond_classes : number of unique bonds

    First perfoms a 1-hot encoding of the bond_type, then passes the
    (*one_hot_encoding, distance) vector to a dense layer. This is the "Edge
    Network" message described by Gilmer, 2017.
    
    """
    
    def __init__(self, units, bond_classes,
                 activation=None,
                 use_bias=True,
                 kernel_initializer='glorot_uniform',
                 bias_initializer='zeros',
                 kernel_regularizer=None,
                 bias_regularizer=None,
                 activity_regularizer=None,
                 kernel_constraint=None,
                 bias_constraint=None,
                 **kwargs):
        
        super(EdgeNetwork, self).__init__(**kwargs)
        self.units = units
        self.bond_classes = bond_classes
        
        self.activation = activations.get(activation)
        self.use_bias = use_bias
        self.kernel_initializer = initializers.get(kernel_initializer)
        self.bias_initializer = initializers.get(bias_initializer)
        self.kernel_regularizer = regularizers.get(kernel_regularizer)
        self.bias_regularizer = regularizers.get(bias_regularizer)
        self.activity_regularizer = regularizers.get(activity_regularizer)
        self.kernel_constraint = constraints.get(kernel_constraint)
        self.bias_constraint = constraints.get(bias_constraint)
        
    def build(self, input_shape):
        
        self.kernel = self.add_weight(shape=(self.bond_classes + 1, self.units**2),
                                      initializer=self.kernel_initializer,
                                      name='kernel',
                                      regularizer=self.kernel_regularizer,
                                      constraint=self.kernel_constraint)
        if self.use_bias:
            self.bias = self.add_weight(shape=(self.units**2,),
                                        initializer=self.bias_initializer,
                                        name='bias',
                                        regularizer=self.bias_regularizer,
                                        constraint=self.bias_constraint)
        else:
            self.bias = None
            
        self.built = True

    def call(self, inputs):
        
        bond_type, distance = inputs
        bond_type_onehot = tf.one_hot(tf.squeeze(bond_type), self.bond_classes)
        stacked_inputs = tf.concat([bond_type_onehot, distance], 1)
        
        output = K.dot(stacked_inputs, self.kernel)
        if self.use_bias:
            output = K.bias_add(output, self.bias, data_format='channels_last')
        if self.activation is not None:
            output = self.activation(output)
            
        output = tf.reshape(output, [-1, self.units, self.units])
        return output
    
    def compute_output_shape(self, input_shape):
        return (input_shape[0][0], self.units, self.units)
    
    def get_config(self):
        config = {
            'units': self.units,
            'bond_classes': self.bond_classes,
            'activation': activations.serialize(self.activation),
            'use_bias': self.use_bias,
            'kernel_initializer': initializers.serialize(self.kernel_initializer),
            'bias_initializer': initializers.serialize(self.bias_initializer),
            'kernel_regularizer': regularizers.serialize(self.kernel_regularizer),
            'bias_regularizer': regularizers.serialize(self.bias_regularizer),
            'activity_regularizer':
                regularizers.serialize(self.activity_regularizer),
            'kernel_constraint': constraints.serialize(self.kernel_constraint),
            'bias_constraint': constraints.serialize(self.bias_constraint)
        }
        base_config = super(EdgeNetwork, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))
