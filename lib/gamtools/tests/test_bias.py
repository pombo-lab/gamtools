import pandas as pd
import numpy as np
from numpy.testing import assert_almost_equal

from gamtools.bias import discretize_feature, observedOverExpected, get_obs_over_exp

def test_discretize_feature():
    discretized = discretize_feature(pd.DataFrame({'name':range(11), 'feature':range(11)}))
    print(discretized.tolist())
    assert discretized.tolist() == ['bin_-inf_to_1.0', 'bin_-inf_to_1.0', 'bin_1.0_to_2.0',
                           'bin_2.0_to_3.0', 'bin_3.0_to_4.0', 'bin_4.0_to_5.0',
                           'bin_5.0_to_6.0', 'bin_6.0_to_7.0', 'bin_7.0_to_8.0',
                           'bin_8.0_to_9.0', 'bin_9.0_to_inf']

def test_discretize_feature_with_nans():
    df = pd.DataFrame({'name':range(12), 'feature':range(12)})
    df.loc[0.0, 'feature'] = np.NaN
    print(df)
    discretized = discretize_feature(df)
    print(discretized.tolist())
    assert discretized.tolist() == [np.NaN, 'bin_-inf_to_2.0', 'bin_-inf_to_2.0',
                                    'bin_2.0_to_3.0', 'bin_3.0_to_4.0', 'bin_4.0_to_5.0',
                                    'bin_5.0_to_6.0', 'bin_6.0_to_7.0', 'bin_7.0_to_8.0',
                                    'bin_8.0_to_9.0', 'bin_9.0_to_10.0', 'bin_10.0_to_inf']

def test_observedOverExpected():
    arr = np.array([[1., 1., 1., 1., 2., 1.],
                    [1., 1., 1., 1., 2., 2.],
                    [1., 1., 1., 1., 1., 1.],
                    [1., 1., 1., 1., 1., 1.],
                    [2., 2., 1., 1., 1., 1.],
                    [1., 2., 1., 1., 1., 1.]])

    o_e_arr = observedOverExpected(arr)
    print(o_e_arr)

    out_arr = np.array([[1,   1,   1,   0.75,1,   1   ],
                        [1,   1,   1,   1,   1.5, 1   ],
                        [1,   1,   1,   1,   1,   0.75],
                        [0.75,1,   1,   1,   1,   1   ],
                        [1,   1.5, 1,   1,   1,   1   ],
                        [1,   1,   0.75,1,   1,   1   ]])
    assert_almost_equal(o_e_arr, out_arr)

def test_get_obs_over_exp():
    arr = np.array([[1., 1., 1., 1., 2., 1.],
                    [1., 1., 1., 1., 2., 2.],
                    [1., 1., 1., 1., 1., 1.],
                    [1., 1., 1., 1., 1., 1.],
                    [2., 2., 1., 1., 1., 1.],
                    [1., 2., 1., 1., 1., 1.]])

    o_e_arr = get_obs_over_exp(arr)
    print(o_e_arr)

    out_arr = np.array([[0,   0,   1,   0.75,1,   1,  ],
                        [0,   0,   0,   1,   1,   1   ],
                        [1,   0,   0,   0,   1,   0.75],
                        [0.75,1,   0,   0,   0,   1   ],
                        [1,   1,   1,   0,   0,   0   ],
                        [1,   1,   0.75,1,   0,   0   ]])
    assert_almost_equal(o_e_arr, out_arr)
