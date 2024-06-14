
%reset -f

import os
import sys
sys.path.insert(0, r"E:\affinewarp-master\affinewarp-master")

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from scipy.io import loadmat
from scipy.io import savemat
from affinewarp import ShiftWarping, PiecewiseWarping

        
os.chdir(r"E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run")
stim_info = loadmat("20201220_Spk_UP_J00005.mat")
data = stim_info['Spk_UP_bin']

 #       data : ndarray, 3d array (trials x times x features) holding signals to be fit.
# model = ShiftWarping(smoothness_reg_scale=20.0)
model = ShiftWarping(maxlag=0.5, smoothness_reg_scale=10.)
# model = ShiftWarping(maxlag=.3, smoothness_reg_scale=10.)
# model = PiecewiseWarping(n_knots=-1, smoothness_reg_scale=30.)

# Fit the model.
model.fit(data, iterations=100)

data_trans = model.transform(data)
loss_hist  = np.asarray(model.loss_hist)
shifts = model.shifts

os.chdir(r"E:\20191019_Michigan_Comp_NSC\Izhikevich_2006_Neural_Computation\20200910_Izhi_run")
savemat(("20201220_Spk_UP_J00005_trans.mat"),{"data_trans":data_trans})
savemat(("20201220_Spk_UP_J00005_loss_hist.mat"),{"loss_hist":loss_hist})
