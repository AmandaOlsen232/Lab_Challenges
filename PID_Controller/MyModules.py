from matplotlib import rcParams
from cycler import cycler
import numpy as np
import csv

class MyPlot:
    def __init__(self):
        # Define your preferred settings
        self.settings = {
            'font.size': 20,
            'font.family': 'Times New Roman',
            'axes.labelsize': 20,
            'axes.titlesize': 16,
            'legend.fontsize': 14,
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            'text.usetex' : True,
            
            'lines.linewidth': 1.5,
            'axes.linewidth': 2,
            'axes.prop_cycle': cycler(color=['black']),
            
            'figure.figsize': (8, 6),
            'savefig.dpi': 600,
            'figure.dpi': 100,
            
            'axes.grid': True,
            'grid.alpha': 1.0,
            'grid.color': 'gray',
            'grid.linestyle': '--',
            'grid.linewidth': 1
        }
        
        rcParams.update(self.settings)
        
    def update(self):
        rcParams.update(self.settings)

def jacobian(f, x, h=1e-3):
    '''

    Parameters
    ----------
    f : function
        name of a function that takes in the state vector, x
        the output has dimension, m.
    x : nx1 array
        n dimensional state vector.

    Returns
    -------
    J : array
        returns mxn Jacobian matrix.

    '''
    
    dim_x = len(x)
    dim_fx = len(f(x))
    J = np.zeros((dim_fx, dim_x)) #initialize Jacobian to the appropriate dimension
    
    #loop through the functions
    for j in range(dim_fx):
        #loop through the input variables
        for i in range(dim_x):
            eps = np.zeros(dim_x).reshape(dim_x,1)
            eps[i][0] = h
            
            fp = f(x+eps)
            fm = f(x-eps)
            J[j][i] = (fp[j][0] - fm[j][0])/(2*h) #central differencing
    return J
        
def multivariable_newtons_method(f, x, tol=1e-10):
    #currently, this function only allows for a symmetric Jacobian (so len(x)=len(f(x)))
    #in the future, I could try and apply a pseudo inverse or something
    #I also need to update the convergence criteria
    '''
    
    Parameters
    ----------
    f : function
        name of a function that takes in the state vector, x
        the output has dimension, m.
    x : nx1 array
        n dimensional state vector.

    Returns
    -------
    x : array
        resulting state vector.
    r : array
        resulting residuals.
    '''
    #do multivariable Newton's method
    r=1000
    iterations = 0
    while r > tol:
        J = jacobian(f, x)
    
        x = x - np.linalg.inv(J) @ f(x)

        r = abs(f(x))
        iterations+=1
    return x, r, iterations

def read_csv(input_csv):
    '''
    
    Parameters
    ----------
    input_csv : string
        .csv file that you want to read into python

    Returns
    -------
    vel : list
        list of csv data with 1 cell per entry in list
    '''
    with open(input_csv, newline='') as file: #flowfield.csv
        reader = csv.reader(file)
        vel = [row for row in reader]
    return vel
