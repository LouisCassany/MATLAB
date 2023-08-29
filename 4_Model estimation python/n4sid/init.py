import pandas as pd
from nfoursid.kalman import Kalman
from nfoursid.nfoursid import NFourSID
from nfoursid.state_space import StateSpace
import matplotlib.pyplot as plt
import numpy as np

FIGSIZE = 8

df = pd.read_csv('data.csv')

nfoursid = NFourSID(
    df,  # the state-space model can summarize inputs and outputs as a dataframe
    output_columns=["psi"],
    input_columns=["delta"],
    num_block_rows=10
)
nfoursid.subspace_identification()

# fig, ax = plt.subplots()
# nfoursid.plot_eigenvalues(ax)
# fig.tight_layout()
# plt.show(block=True)

ORDER_OF_MODEL_TO_FIT = 4
state_space_identified, covariance_matrix = nfoursid.system_identification(
    rank=ORDER_OF_MODEL_TO_FIT
)

print(state_space_identified.a)