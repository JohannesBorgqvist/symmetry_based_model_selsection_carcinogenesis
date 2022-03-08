import matplotlib.pyplot as plt
import pandas as pd
from hampel import hampel
import numpy as np
# Define the y-value
y = np.array([1, 2, 1 , 1 , 1, 2, 13, 2, 1, 2, 15, 1, 2])/10
x = np.arange(len(y))

ts = pd.Series(list(y))

# Just outlier detection
outlier_indices = hampel(ts, window_size=5, n=3)
print("Outlier Indices: ", outlier_indices)

new_x = np.delete(x, outlier_indices)
new_y = np.delete(y, outlier_indices)

# Outlier Imputation with rolling median
#ts_imputation = hampel(ts, window_size=5, n=3, imputation=True)


fig, axes = plt.subplots(1,2,figsize=(15,5))
axes[0].plot(x,y,label='Original data')
axes[0].legend()
axes[1].plot(new_x,new_y,label='Removed outliers')
axes[1].legend()
fig.add_subplot(111, frameon=False)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Age, $t$")
plt.ylabel("Incidence, $R(t)$")
# displaying the title
plt.title("Fit of candidate models to three data sets",fontsize=20, fontweight='bold')
plt.show()
