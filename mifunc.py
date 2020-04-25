import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from filterpy.kalman import KalmanFilter


def tracker1(x_initial, R_std, Q_std):
    tracker = KalmanFilter(dim_x=1, dim_z=1)

    tracker.F = np.array([1])
    tracker.u = 0.
    tracker.H = np.array([1])

    tracker.R = R_std
    tracker.Q = Q_std
    tracker.x = np.array([x_initial]).T
    tracker.P = 50

    return tracker


def rssi_filter(data, R_std, Q_std):

    x_initial = data[0]

    # run filter
    robot_tracker = tracker1(x_initial, R_std, Q_std)
    mu, cov, _, _ = robot_tracker.batch_filter(data)

    res_nan = np.isnan(mu)
    mu = mu[~res_nan]
    mu = mu.tolist()
    
    return mu
