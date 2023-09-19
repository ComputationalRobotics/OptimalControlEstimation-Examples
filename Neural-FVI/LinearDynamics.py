import numpy as np
import torch as th

class LinearSystem:
    def __init__(self, dt = 0.01):
        self.dt = dt
    
    def dynamics(self, x, u):
        # Double Integrator
        self.fx = np.array([x[1], 0])
        self.gx = np.array([0, 1])
        xdot = self.fx + self.gx * u
        return self.fx, self.gx, xdot
    
    def dynamics_torch(self, x_batch, device):
        x0, x1 = x_batch[:, 0], x_batch[:, 1]
        zero_tensor = th.zeros_like(x1)
        fx = th.stack((x1, zero_tensor), dim=1)
        gx = th.tensor([0, 1], dtype=th.float32, device=device, requires_grad=True).repeat(x_batch.size(0), 1)
        return fx, gx
    
