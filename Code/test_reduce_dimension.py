import matplotlib.pyplot as plt
import numpy as np
import clean_up_output

x_dense = np.linspace(-2,2,500,endpoint=True)
y_dense = x_dense**2+2
x_sparse = np.linspace(-2,2,100,endpoint=True)
y_sparse = clean_up_output.reduce_density(x_dense,y_dense,x_sparse)


fig, ax = plt.subplots()
ax.plot(x_dense, y_dense,label="Dense function, $y=x^2+2$")
ax.plot(x_sparse, y_sparse,label="Sparse function, $y=x^2+2$")
ax.set(xlabel='Variable $x$', ylabel='Function value, $y(x)$',
       title='Plot of densities')
ax.legend()
ax.grid()
plt.show()
