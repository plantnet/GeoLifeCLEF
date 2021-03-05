import numbers

import numpy as np


def top_k_error_rate(y_true, y_score, k):
    r"""Computes the top_k error rate for a given k.

    Parameters
    ----------
    y_true: 1d array, [n_samples]
        True labels.
    y_score: 2d array, [n_samples, n_classes]
        Scores for each label.
    k: int
        Value of k to use, should range from 1 to n_classes (default: 30).

    Returns
    -------
    float:
        Error rate value.

    Notes
    -----
    Complexity: :math:`O( n_\text{samples} \times n_\text{classes} )`.
    """
    # y_true = validate_labels(y_true)
    # y_score = validate_scores(y_score)

    n_classes = y_score.shape[1]

    if not (isinstance(k, numbers.Integral) and 0 < k <= n_classes):
        raise ValueError(
            "k should be an integer ranging from 1 to n_classes, given {}".format(k)
        )

    # Compute predictions
    y_pred = np.argpartition(y_score, n_classes - k, axis=1)[:, -k:]

    # Compute error rate
    pointwise_accuracy = np.sum(y_pred == y_true[:, None], axis=1)
    return 1 - np.mean(pointwise_accuracy)


def top_30_error_rate(y_true, y_score):
    r"""Computes the top_k error rate for a given k=30.

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
