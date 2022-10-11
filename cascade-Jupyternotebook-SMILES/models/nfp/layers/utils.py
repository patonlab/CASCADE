""" Unused for the moment, I believe. """

import tensorflow as tf

def get_shape(tensor):
    """Returns the tensor's shape.
    Each shape element is either:
    - an `int`, when static shape values are available, or
    - a `tf.Tensor`, when the shape is dynamic.
    Args:
    tensor: A `tf.Tensor` to get the shape of.
    Returns:
    The `list` which contains the tensor's shape.
    """

    shape_list = tensor.shape.as_list()
    if all(s is not None for s in shape_list):
        return shape_list
    
    shape_tensor = tf.shape(tensor)
    return [shape_tensor[i] if s is None else s for i, s in
            enumerate(shape_list)]


def repeat(tensor, repeats, axis=0):
    """Repeats a `tf.Tensor`'s elements along an axis by custom amounts.
    Equivalent to Numpy's `np.repeat`.
    `tensor and `repeats` must have the same numbers of elements along `axis`.
    Args:
    tensor: A `tf.Tensor` to repeat.
    repeats: A 1D sequence of the number of repeats per element.
    axis: An axis to repeat along. Defaults to 0.
    name: (string, optional) A name for the operation.
    Returns:
    The `tf.Tensor` with repeated values.
    """

    cumsum = tf.cumsum(repeats)
    range_ = tf.range(cumsum[-1])

    indicator_matrix = tf.cast(tf.expand_dims(range_, 1) >= cumsum, tf.int32)
    indices = tf.reduce_sum(indicator_matrix, reduction_indices=1)

    shifted_tensor = _axis_to_inside(tensor, axis)
    repeated_shifted_tensor = tf.gather(shifted_tensor, indices)
    repeated_tensor = _inside_to_axis(repeated_shifted_tensor, axis)

    shape = tensor.shape.as_list()
    shape[axis] = None
    repeated_tensor.set_shape(shape)

    return repeated_tensor


def _axis_to_inside(tensor, axis):
    """Shifts a given axis of a tensor to be the innermost axis.
    Args:
        tensor: A `tf.Tensor` to shift.
        axis: An `int` or `tf.Tensor` that indicates which axis to shift.
    Returns:
        The shifted tensor.
    """

    axis = tf.convert_to_tensor(axis)
    rank = tf.rank(tensor)

    range0 = tf.range(0, limit=axis)
    range1 = tf.range(tf.add(axis, 1), limit=rank)
    perm = tf.concat([[axis], range0, range1], 0)

    return tf.transpose(tensor, perm=perm)


def _inside_to_axis(tensor, axis):
    """Shifts the innermost axis of a tensor to some other axis.
    Args:
        tensor: A `tf.Tensor` to shift.
        axis: An `int` or `tf.Tensor` that indicates which axis to shift.
    Returns:
        The shifted tensor.
    """

    axis = tf.convert_to_tensor(axis)
    rank = tf.rank(tensor)

    range0 = tf.range(1, limit=axis + 1)
    range1 = tf.range(tf.add(axis, 1), limit=rank)
    perm = tf.concat([range0, [0], range1], 0)

    return tf.transpose(tensor, perm=perm)
