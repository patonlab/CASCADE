from keras.models import Model

class GraphModel(Model):
    """ This is a simple modification of the Keras `Model` class to avoid
    checking each input for a consistent batch_size dimension. Should work as
    of keras-team/keras#11548.

    """

    def _standardize_user_data(self, *args, **kwargs):
        kwargs['check_array_lengths'] = False
        return super(GraphModel, self)._standardize_user_data(*args, **kwargs)
