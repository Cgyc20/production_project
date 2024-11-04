import numpy as np
from tqdm import tqdm
import os
import json

class PDE:
    def __init__(self, domain_length, PDE_points, total_time, timestep, production_rate, degradation_rate, diffusion_rate, PDE_initial):
        self.L = domain_length
        self.PDE_points = PDE_points
        self.production_rate = production_rate
        self.deltax = self.L / self.PDE_points
        self.total_time = total_time
        self.PDE_initial_conditions = PDE_initial
        self.timestep = timestep
        self.diffusion_rate = diffusion_rate
        self.degradation_rate = degradation_rate
        
        self.PDE_X = np.linspace(0, self.L - self.deltax, self.PDE_points)
        self.steady_state = production_rate / degradation_rate
        self.DX_NEW = self.create_finite_difference()
        self.time_vector = np.arange(0, total_time, timestep)
        self.Crank_matrix, self.M1_inverse = self.create_crank_nicholson()

        self.PDE_grid = np.zeros((self.PDE_points, len(self.time_vector)))
        self.PDE_grid[:, 0] = self.PDE_initial_conditions

        print("Successfully initialized the hybrid model")

    def create_crank_nicholson(self):
        H = self.create_finite_difference()
        M1 = np.identity(H.shape[0]) * (1 + 0.5 * self.timestep * self.degradation_rate) - 0.5 * (self.timestep * self.diffusion_rate / self.deltax**2) * H
        M2 = np.identity(H.shape[0]) * (1 - 0.5 * self.timestep * self.degradation_rate) + 0.5 * (self.timestep * self.diffusion_rate / self.deltax**2) * H
        M1_inverse = np.linalg.inv(M1)
        Crank_matrix = M1_inverse @ M2
        return Crank_matrix, M1_inverse

    def create_finite_difference(self):
        self.DX = np.zeros((self.PDE_points, self.PDE_points), dtype=int)
        self.DX[0, 0], self.DX[-1, -1] = -1, -1
        self.DX[0, 1], self.DX[-1, -2] = 1, 1
        for i in range(1, self.DX.shape[0] - 1):
            self.DX[i, i] = -2
            self.DX[i, (i + 1)] = 1
            self.DX[i, (i - 1)] = 1
        return self.DX

    def crank_nicholson(self, old_vector):
        return self.Crank_matrix @ old_vector + self.M1_inverse@ (self.production_rate * self.timestep*np.ones(self.PDE_points))

    def run_simulation(self):
        for i in range(len(self.time_vector) - 1):
            self.PDE_grid[:, i + 1] = self.crank_nicholson(self.PDE_grid[:, i])
        print("Simulation completed")
        return self.PDE_grid

    def save_simulation_data(self, PDE_grid, datadirectory='data'):
        if not os.path.exists(datadirectory):
            os.makedirs(datadirectory)
        params = {
            'domain_length': self.L,
            'PDE_points': self.PDE_points,
            'total_time': self.total_time,
            'timestep': self.timestep,
            'production_rate': self.production_rate,
            'degradation_rate': self.degradation_rate,
            'diffusion_rate': self.diffusion_rate,
        }
        np.savez(os.path.join(datadirectory, "PDE_data.npz"), PDE_grid=PDE_grid, PDE_X=self.PDE_X, time_vector=self.time_vector, parameters=params)
        print("Data saved successfully")

