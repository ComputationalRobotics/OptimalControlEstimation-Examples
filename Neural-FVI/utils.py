import numpy as np
import scipy.linalg
from scipy.integrate import odeint
from control import ctrb
import matplotlib.pyplot as plt
import torch as th
from torch.func import jacrev, hessian, vmap

def sample_state(x0max, x1max):
    x0 = (-1 + 2 * np.random.rand()) * x0max
    x1 = (-1 + 2 * np.random.rand()) * x1max
    
    return(np.array([x0, x1]))

def prepare_data(sample_number, x0max, x1max, vector_length=2):
    x_train = np.empty((sample_number,vector_length))
    for i in range(sample_number):
        new_element = sample_state(x0max=x0max, x1max=x1max)
        x_train[i] = new_element
        
    return x_train
    

if __name__ == "__main__":
    pass