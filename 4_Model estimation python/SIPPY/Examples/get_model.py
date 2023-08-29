from __future__ import division

from past.utils import old_div

try:
    from sippy import *
except ImportError:
    import sys
    import os

    sys.path.append(os.pardir)
    from sippy import *

import numpy as np
from sippy import functionset as fset
from sippy import functionsetSIM as fsetSIM
import matplotlib.pyplot as plt

# import os
# os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

data = np.genfromtxt('data.csv', delimiter=',')
data = data[1:]

U = data[:,4]
y_tot = data[:,3]
Time = data[:,0]
U = [U]
y_tot = [y_tot]

plt.close("all")
plt.figure(0)
plt.plot(Time, U[0])
plt.ylabel("input")
plt.grid()
plt.xlabel("Time")
#
plt.figure(1)
plt.plot(Time, y_tot[0])
plt.ylabel("y_tot")
plt.grid()
plt.xlabel("Time")
plt.title("Ytot")

METHOD = ['N4SID', 'CVA', 'MOESP', 'PARSIM-S', 'PARSIM-P', 'PARSIM-K']

#System identification
METHOD = ['N4SID', 'CVA', 'MOESP', 'PARSIM-S', 'PARSIM-P', 'PARSIM-K']
lege = ['System']
for i in range(len(METHOD)):
    method = METHOD[i]
    sys_id = system_identification(y_tot, U, method, SS_fixed_order = 2 )
    print(sys_id.x0)
    xid, yid = fsetSIM.SS_lsim_process_form(sys_id.A, sys_id.B, sys_id.C, sys_id.D, U, sys_id.x0)
    #
    plt.plot(Time, yid[0])
    plt.show()
    lege.append(method) 
plt.legend(lege) 


