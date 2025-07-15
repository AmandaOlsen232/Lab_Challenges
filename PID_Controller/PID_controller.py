import math
import matplotlib.pyplot as plt
import numpy as np
import MyModules as my
my.MyPlot()

class PID:
    
    def __init__(self, solver, kp=-1, ki=-5, kd=0, p_final=0.5):
        self.V = 520.0 #ft/s
        self.rho = 0.00187258 #slug/ft^3
        self.altitude = 8000.0 #ft
        
        self.Sw = 300.0 #ft^2
        self.bw = 30.0 #ft
        self.weight = 20500.0 #lbf
        self.Ixx = 9496.0 #slug-ft^2

        self.Cl_delta_a = -0.1172
        self.Cl_pbar = -0.375
        self.sigma_r = -self.rho*self.Sw*(self.bw**2)*self.V*self.Cl_pbar/(4*self.Ixx)
        
        self.p_final = p_final
        self.t_final = 3.0
        
        self.delta_t = 0.01
        
        #starting values at time t=0
        self.t = 0.0 
        self.p = 0.0
        self.error = self.p_final - self.p
        self.delta_a = 0.0
        self.int_error = 0.0
        
        self.kp = kp
        self.ki = ki
        self.kd = kd
        
        self.solver = solver
        
        self.track_error = []
        self.track_time = []
        self.track_p = []
        self.track_delta_a = []      


        
    def P(self):
        return self.kp*self.error
    
    def I(self):
        count = 0
        self.int_error = 0.0
        for _ in range(len(self.track_error)-1):
            f1 = self.track_error[count]
            f2 = self.track_error[count+1]
            self.int_error += self.delta_t*(f1+f2)/2
            count += 1
        
        return self.int_error*self.ki
    
    def D(self):
        return self.kd*((self.track_error[-1] - self.track_error[-2])/self.delta_t)
        
    def update_error(self):
        self.error = self.p_final - self.p
    
    def update_time(self):
        self.t += self.delta_t

    def update_p(self):
        if self.solver=="analytic":
            self.p = (-2*self.V/self.bw)*(self.Cl_delta_a*self.delta_a/self.Cl_pbar)*(1 - math.exp(-self.sigma_r*self.t))
        if self.solver=="numerical":
            # dp_dt = (-2*self.V/self.bw)*(self.Cl_delta_a*self.delta_a/self.Cl_pbar)*(self.sigma_r*math.exp(-self.sigma_r*self.t))
            # self.p += dp_dt*self.delta_t
            pbar = self.p*self.bw/(2*self.V)
            dp_dt = (self.rho*self.Sw*self.bw*(self.V**2)/(2*self.Ixx))*(self.Cl_pbar*pbar + self.Cl_delta_a*self.delta_a)
            self.p += dp_dt*self.delta_t
            
    def update_delta_a(self):
        self.delta_a = self.P() + self.I()
    
    def static_delta_a(self):
        self.delta_a = -self.p_final*self.bw*self.Cl_pbar/(2*self.V*self.Cl_delta_a)
    
    def track(self):
        self.track_error.append(self.error)
        self.track_time.append(self.t)
        self.track_p.append(self.p)
        self.track_delta_a.append(self.delta_a)
        
if __name__ == "__main__":
    n = PID("numerical", kp=-0.5, ki=-4)
    n.track()
    
    for _ in range(int(n.t_final/n.delta_t)):
        n.update_time()
        n.update_delta_a()
        n.update_p()
        n.update_error()
        n.track()
        
    error = np.array(n.track_error)
    time = np.array(n.track_time)
    p = np.array(n.track_p)
    delta_a = np.array(n.track_delta_a)
    
    plt.plot(time, p, color="blue", label="PID Controller")
    plt.xlabel("Time [sec]")
    plt.ylabel("p")
    plt.legend()
    
    # a = PID("analytic")
    # a.static_delta_a()
    # a.track()
    
    # for _ in range(int(a.t_final/a.delta_t)):
    #     a.update_time()
    #     # n.update_delta_a()
    #     a.update_p()
    #     a.update_error()
    #     a.track()
        
    # error = np.array(a.track_error)
    # time = np.array(a.track_time)
    # p = np.array(a.track_p)
    # delta_a = np.array(a.track_delta_a)
    
    # plt.plot(time, p, color="red", label="analytic")
    