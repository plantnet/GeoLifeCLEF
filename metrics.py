import numbers
import warnings

import numpy as np


def generic_validation(
    arr, name, ndim=None, dtype=None, allow_nan=True, allow_inf=True
):
    """Generic array validation

    Check if the given input satisfies the given constraints.

    Parameters
    ----------
    arr : array-like
        Array to validate.
    name : string
        Name of the array to use in exceptions.
    ndim : integer
        Expected array number of dimensions.
    dtype : data type
        Expected array data type.
    allow_nan : boolean
        If False, checks there is no NaN in the array.
    allow_inf : boolean
        If False, check there is no inf in the array.

    Returns
    -------
    arr : array
        Validated array.

    Raises
    ------
    ValueError
        If the array does not satisfy the constraints.
    """
    shape = arr.shape

    if ndim is not None and arr.ndim != ndim:
        raise ValueError(
            "{} should be a {}d array, given array of shape {}".format(
                name, ndim, shape
            )
        )

    if dtype is not None and not np.issubdtype(arr.dtype, dtype):
        raise ValueError(
            "{} should be of type {} (or derived), given {}".format(
                name, dtype, arr.dtype
            )
        )

    if not allow_nan and np.isnan(arr).any():
        raise ValueError("{} should only contain non-NaN values".format(name))

    if not allow_inf and not np.isfinite(arr).all():
        raise ValueError("{} should only contain finite values".format(name))

    return arr


def validate_labels(labels):
    """Validate labels array

    Check if the given input is a 1d array with finite non-NaN integers values.

    Parameters
    ----------
    labels : array-like
        Candidate label array to validate.

    Returns
    -------
    labels : array
        Validated array.

    Raises
    ------
    ValueError
        If the array does not satisfy the constraints.
    """
    labels = np.asarray(labels)
    return generic_validation(
        labels, "Labels", ndim=1, dtype=np.integer, allow_inf=False, allow_nan=False
    )


def validate_scores(scores):
    """Validate scores array

    Check if the given input is a 2d array with non-NaN numerical values.

    Parameters
    ----------
    scores : array-like
        Candidate label array to validate.

    Returns
    -------
    scores : array
        Validated array.

    Raises
    ------
    ValueError
        If the array does not satisfy the constraints.
    """
    scores = np.asarray(scores)
    return generic_validation(
        scores, "Scores", ndim=2, dtype=np.number, allow_inf=True, allow_nan=False
    )


def validate_top_k_sets(s_pred):
    """Validate top-k sets

    Check if the given input is an array representing sets.

    Parameters
    ----------
    s_pred : array
        Candidate sets to validate.

    Returns
    -------
    s_pred : array
        Validated sets.

    Raises
    ------
    ValueError
        If the array does not satisfy the constraints.
    """
    s_pred = np.asarray(s_pred)
    return generic_validation(
        s_pred, "Top-k sets", ndim=2, dtype=np.integer, allow_inf=False, allow_nan=False
    )


def predict_top_k_set(y_score, k, disable_warning=False):
    r"""Predicts the top-k sets from scores for a given k.

    Parameters
    ----------
    y_score: 2d array, [n_samples, n_classes]
        Scores for each sample and label.
    k: int
        Value of k to use, should range from 1 to n_classes.
    disable_warning: bool
        Disables the warning trigger if y_score looks like it contains top-k sets rather than scores

    Returns
    -------
    2d array, [n_samples, k]:
        Predicted top-k sets for each sample.

    Notes
    -----
    Complexity: :math:`O( n_\text{samples} \times n_\text{classes} )`.
    """
    if not disable_warning:
        try:
            y_score = validate_top_k_sets(y_score)

            if np.issubdtype(y_score.dtype, np.integer) and y_score.shape[1] == k:
                warnings.warn(
                    "y_score is an integer array with already {} columns".format(k),
                    UserWarning,
                )
        except ValueError:
            pass

    y_score = validate_scores(y_score)

    n_classes = y_score.shape[1]

    if not (isinstance(k, numbers.Integral) and 0 < k <= n_classes):
        raise ValueError(
            "k should be an integer ranging from 1 to n_classes, given {}".format(k)
        )

    n_classes = y_score.shape[1]
    s_pred = np.argpartition(y_score, n_classes - k, axis=1)[:, -k:]

    # Return a copy instead of a view to save memory
    return s_pred.copy()


def predict_top_30_set(y_score):
    r"""Predicts the top-30 sets from scores.

    Parameters
    ----------
    y_score: 2d array, [n_samples, n_classes]
        Scores for each sample and label.

    Returns
    -------
    2d array, [n_samples, 30]:
        Predicted top-30 sets for each sample.

    Notes
    -----
    Complexity: :math:`O( n_\text{samples} \times n_\text{classes} )`.
    """
    return predict_top_k_set(y_score, k=30)


def top_k_error_rate_from_sets(y_true, s_pred):
    r"""Computes the top-k error rate from predicted sets.

    Parameters
    ----------
    y_true: 1d array, [n_samples]
        True labels.
    s_pred: 2d array, [n_samples, k]
        Previously computed top-k sets for each sample.

    Returns
    -------
    float:
        Error rate value.
    """
    y_true = validate_labels(y_true)
    s_pred = validate_top_k_sets(s_pred)

    pointwise_accuracy = np.sum(s_pred == y_true[:, None], axis=1)
    return 1 - np.mean(pointwise_accuracy)


def top_k_error_rate(y_true, y_score, k, disable_warning=False):
    r"""Computes the top-k error rate for a given k.

    Parameters
    ----------
    y_true: 1d array, [n_samples]
        True labels.
    y_score: 2d array, [n_samples, n_classes]
        Scores for each label.
    k: int
        Value of k to use, should range from 1 to n_classes.
    disable_warning: bool
        Disables the warning trigger if y_score looks like it contains top-k sets rather than scores

    Returns
    -------
    float:
        Error rate value.

    Notes
    -----
    Complexity: :math:`O( n_\text{samples} \times n_\text{classes} )`.
    """
    s_pred = predict_top_k_set(y_score, k, disable_warning=disable_warning)
    return top_k_error_rate_from_sets(y_true, s_pred)


def top_30_error_rate(y_true, y_score):
    r"""Computes the top-30 error rate.

    Parameters
    ----------
    y_true: 1d array, [n_samples]
        True labels.
    y_score: 2d array, [n_samples, n_classes]
        Scores for each label.

    Returns
    -------
    float:
        Top-30 error rate value.

    Notes
    -----
    Complexity: :math:`O( n_\text{samples} \times n_\text{classes} )`.
    """
    return top_k_error_rate(y_true, y_score, k=30)
