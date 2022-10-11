"""
These classes wrap the GRUCell and LSTMCell keras layers, allowing weight-tying
between RNNs for each atom. Essentially, Keras expects inputs to an RNN to be
of the shape (batch_size, sequence_length, features), while our inputs are of
shape (num_atoms_per_batch, features). 

Peter, 3/20/2018

I should revisit whether these are needed post graph_nets update.

Peter, 11/9/2018.

"""

from keras.layers import GRUCell, LSTMCell

class GRUStep(GRUCell): 
    def build(self, input_shape):
        GRUCell.build(self, input_shape[0])

    def call(self, inputs):
        """ Inputs should be [message, previous_state], returns [next_state]
        """
        return GRUCell.call(self, inputs[0], [inputs[1]])[0]

    def compute_output_shape(self, input_shape):
        return input_shape[0]


class LSTMStep(LSTMCell): 
    def build(self, input_shape):
        LSTMCell.build(self, input_shape[0])

    def call(self, inputs):
        """ Inputs should be [message, previous_state, previous_memory], 
        returns [next_state, next_memory]
        """
        outputs = LSTMCell.call(self, inputs[0], [inputs[1], inputs[2]])
        return [outputs[0], outputs[1][1]]

    def compute_mask(self, inputs, mask):
        return [None] * 2
        
    def compute_output_shape(self, input_shape):
        return [input_shape[0], input_shape[0]]
