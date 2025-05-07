import numpy as np
import os

class PDE:
    def __init__(self, domain_length, PDE_points, total_time, timestep,
                 production_rate, degradation_rate, diffusion_rate, PDE_initial):
        
        self.L = domain_length
        self.PDE_points = PDE_points
        self.deltax = self.L / self.PDE_points
        self.total_time = total_time
        self.timestep = timestep
        self.production_rate = production_rate
        self.degradation_rate = degradation_rate
        self.diffusion_rate = diffusion_rate

        self.PDE_X = np.linspace(0, self.L, self.PDE_points)
        self.time_vector = np.arange(0, total_time + timestep, timestep)
        self.PDE_grid = np.zeros((self.PDE_points, len(self.time_vector)))
        self.PDE_grid[:, 0] = PDE_initial

        self.H = self.create_finite_difference()

        print("Successfully initialized the PDE model using RK4")

    def create_finite_difference(self):
        H = np.zeros((self.PDE_points, self.PDE_points))
        H[0, 0], H[-1, -1] = -1, -1
        H[0, 1], H[-1, -2] = 1, 1
        for i in range(1, self.PDE_points - 1):
            H[i, i - 1] = 1
            H[i, i] = -2
            H[i, i + 1] = 1
        return H

    def RHS_spatial_deriv(self, u):
        production_vector = np.zeros_like(u)
        production_vector[0] = self.production_rate / self.deltax  # Neumann BC
        degradation_term = self.degradation_rate * u
        diffusion_term = (self.diffusion_rate / self.deltax**2) * self.H @ u
        return diffusion_term + production_vector - degradation_term

    def runge_kutta_4(self, u):
        k1 = self.RHS_spatial_deriv(u)
        k2 = self.RHS_spatial_deriv(u + 0.5 * self.timestep * k1)
        k3 = self.RHS_spatial_deriv(u + 0.5 * self.timestep * k2)
        k4 = self.RHS_spatial_deriv(u + self.timestep * k3)
        return u + (self.timestep / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

    def run_simulation(self):
        for i in range(len(self.time_vector) - 1):
            self.PDE_grid[:, i + 1] = self.runge_kutta_4(self.PDE_grid[:, i])
        print("Simulation completed")
        return self.PDE_grid

    def save_simulation_data(self,PDE_grid,datadirectory='data'):
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
        np.savez(os.path.join(datadirectory, "PDE_data.npz"),
                 PDE_grid=PDE_grid,
                 PDE_X=self.PDE_X,
                 time_vector=self.time_vector,
                 parameters=params)
        print("Data saved successfully")

