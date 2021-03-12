import numpy as np
import pytest

from metrics import validate_top_k_sets, predict_top_k_set, top_k_error_rate_from_sets, top_k_error_rate


def assert_same_top_k_sets(result, expected):
    result = validate_top_k_sets(result)
    expected = validate_top_k_sets(expected)

    result = np.sort(result, axis=1)
    expected = np.sort(expected, axis=1)

    assert result.shape == expected.shape
    assert np.all(result == expected)


def test_predict_top_k_set():
    y_score = [
        [1, 2, 3, 4, 5],
        [5, 4, 3, 2, 1],
        [2, 1, 3, 5, 4],
    ]

    k = 1
    result = predict_top_k_set(y_score, k=k)
    expected = np.asarray([
        [4],
        [0],
        [3],
    ])
    assert_same_top_k_sets(result, expected)

    k = 3
    result = predict_top_k_set(y_score, k=k)
    expected = np.asarray([
        [4, 3, 2],
        [0, 1, 2],
        [3, 4, 2],
    ])
    assert_same_top_k_sets(result, expected)

    y_score = np.array(y_score, dtype=np.float32)
    result = predict_top_k_set(y_score, k=k)
    assert_same_top_k_sets(result, expected)

    with pytest.warns(UserWarning):
        predict_top_k_set(result, k)


def test_top_k_error_rate_from_sets():
    y_true = [0, 1, 2]
    y_score = np.asarray([
        [1, 2, 3, 4, 5],
        [5, 4, 3, 2, 1],
        [2, 1, 3, 5, 4],
    ], dtype=np.float32)
    s_pred = predict_top_k_set(y_score, k=3)

    with pytest.raises(ValueError):
        top_k_error_rate_from_sets(y_true, y_score)

    result = top_k_error_rate_from_sets(y_true, s_pred)
    expected = (1 + 0 + 0) / 3
    assert result == pytest.approx(expected)


def test_top_k_error_rate():
    k = 3
    y_true = [0, 1, 2]
    y_score = np.asarray([
        [1, 2, 3, 4, 5],
        [5, 4, 3, 2, 1],
        [2, 1, 3, 5, 4],
    ], dtype=np.float32)

    s_pred = predict_top_k_set(y_score, k)
    with pytest.warns(UserWarning):
        top_k_error_rate(y_true, s_pred, k)

    result = top_k_error_rate(y_true, y_score, k)
    expected = (1 + 0 + 0) / 3
    assert result == pytest.approx(expected)
