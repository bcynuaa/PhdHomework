import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import copy

rk4_coeff = np.array([
    1/4, 1/3, 1/2, 1
])

###################################################################################################

class Case:
    
    def __init__(self, Ma_inf, alpha):
        self.Ma_inf = Ma_inf
        self.alpha = alpha * np.pi / 180
        self.T_inf = 1.0
        self.C_inf = 1.0
        self.P_inf = 1.0
        self.rho_inf = 1.4
        self.R = 1 / 1.4
        self.gamma = 1.4
        self.Cv = self.R / (self.gamma - 1.0)
        self.calculate_conserve_inf()
        self.k2 = 0.8
        self.k4 = 3e-3
        self.CFL = 2.0
        self.error = 1e-5
        self.STEP = 10
        self.data_path = ".//data//"
        self.image_path = ".//image//"
        pass
    
    def calculate_conserve_inf(self):
        w_inf = np.zeros(4)
        w_inf[0] = self.rho_inf
        U_inf = self.Ma_inf * self.C_inf * np.cos(self.alpha)
        V_inf = self.Ma_inf * self.C_inf * np.sin(self.alpha)
        w_inf[1] = self.rho_inf * U_inf
        w_inf[2] = self.rho_inf * V_inf
        E_inf = self.P_inf / (self.gamma - 1) / self.rho_inf + (U_inf**2 + V_inf**2) / 2
        w_inf[3] = self.rho_inf * E_inf
        self.w_inf = w_inf
        self.w_inf_all = self.conserve_to_all(w_inf)
        pass
    
    def conserve_to_all(self, w):
        rho = w[0]
        U = w[1] / rho
        V = w[2] / rho
        E = w[3] / rho
        P = (E - (U**2 + V**2) / 2) * rho * (self.gamma - 1)
        H = E + P / rho
        C = np.sqrt( self.gamma * P / rho )
        return np.array([rho, U, V, E, P, H, C])
        pass
    
    def calculate_flux(self, w_edge, dx, dy):
        rho, U, V, E, P, H, C = self.conserve_to_all(w_edge)
        Q = np.zeros(4)
        Z = U*dy - V*dx
        Q[0] = Z * rho
        Q[1] = Z * rho * U + P * dy
        Q[2] = Z * rho * V - P * dx
        Q[3] = Z * rho * H
        return Q
        pass
    pass

###################################################################################################

class Grid:
    
    def __init__(self, grid_file):
        nnodes, nedges, ncells = np.loadtxt(grid_file, max_rows=1, dtype=int)
        xy = np.loadtxt(grid_file, skiprows=1, max_rows=nnodes, dtype=float)
        iedge = np.loadtxt(grid_file, skiprows=1+nnodes, max_rows=nedges, dtype=int) - 1
        icell = np.loadtxt(grid_file, skiprows=1+nnodes+nedges, max_rows=ncells, dtype=int) - 1
        vol = np.loadtxt(grid_file, skiprows=1+nnodes+nedges+ncells, max_rows=ncells, dtype=float)
        self.nnodes = nnodes
        self.nedges = nedges
        self.ncells = ncells
        self.xy = xy
        self.iedge = iedge
        self.icell = icell
        self.vol = vol
        self.calculate_dx_dy()
        pass
    
    def plot_grid(self, image_name):
        triangle = tri.Triangulation(self.xy.T[0], self.xy.T[1], self.icell)
        plt.figure(figsize=(12, 8), facecolor="white", dpi=500)
        plt.gca().set_aspect(1)
        plt.triplot(triangle, lw=0.1, color="r")
        plt.xlim(-2, 3)
        plt.ylim(-2.5, 2.5)
        plt.savefig(image_name, bbox_inches="tight")
        plt.show()
        pass
    
    def calculate_dx_dy(self):
        dx = np.zeros(self.nedges)
        dy = np.zeros(self.nedges)
        ds = np.zeros(self.nedges)
        for edge in range(self.nedges):
            i, j, k, p = self.iedge[edge]
            xj, yj = self.xy[j]
            xi, yi = self.xy[i]
            dx[edge] = xj - xi
            dy[edge] = yj - yi
            ds[edge] = np.sqrt( dx[edge] ** 2  + dy[edge]**2 )
            pass
        self.dx = dx
        self.dy = dy
        self.ds = ds
        pass
    
    pass

###################################################################################################

class Solver:
    
    def __init__(self, case, grid):
        self.case = case
        self.grid = grid
        pass
    
    def initialize_field(self):
        W_cur = np.zeros([self.grid.ncells, 4])
        for j in range(self.grid.ncells):
            W_cur[j] = self.case.w_inf
            pass
        self.W_cur = W_cur
        pass
    
    def conserve_to_all(self):
        W_cur_all = np.zeros([self.grid.ncells, 7])
        for j in range(self.grid.ncells):
            W_cur_all[j] = self.case.conserve_to_all(self.W_cur[j])
            pass
        self.W_cur_all = W_cur_all
        pass
    
    def calculate_flux(self, rho, U, V, P, H, dx, dy):
        Q = np.zeros(4)
        Z = U*dy - V*dx
        Q[0] = Z * rho
        Q[1] = Z * rho * U + P * dy
        Q[2] = Z * rho * V - P * dx
        Q[3] = Z * rho * H
        return Q
        pass
    
    def calculate_laplace(self):
        W_cur_laplace = np.zeros([self.grid.ncells, 4])
        for edge in range(self.grid.nedges):
            i, j, k, p = self.grid.iedge[edge]
            if p > -1:
                w_k = self.W_cur[k]
                w_p = self.W_cur[p]
                value = w_p - w_k
                W_cur_laplace[k] = W_cur_laplace[k] + value
                W_cur_laplace[p] = W_cur_laplace[p] - value
                pass
            pass
        self.W_cur_laplace = W_cur_laplace
        pass
    
    def each_step(self):
        self.conserve_to_all()
        self.calculate_laplace()
        Q = np.zeros([self.grid.ncells, 4])
        D2 = np.zeros([self.grid.ncells, 4])
        D4 = np.zeros([self.grid.ncells, 4])
        t = np.zeros(self.grid.ncells)
        for edge in range(self.grid.nedges):
            i, j, k, p = self.grid.iedge[edge]
            dx = self.grid.dx[edge]
            dy = self.grid.dy[edge]
            ds = self.grid.ds[edge]
            if p == -2:
                w_k = self.W_cur[k]
                p = self.W_cur_all[k][4]
                flux = np.array([
                    0, p * dy, -p * dx, 0
                ])
                Q[k] = Q[k] + flux
                pass
            elif p == -3:
                w_k = self.W_cur[k]
                w_k_all = self.W_cur_all[k]
                w_inf = self.case.w_inf
                w_inf_all = self.case.w_inf_all
                
                C_k = w_k_all[6]
                rho_k = w_k_all[0]
                P_k = w_k_all[4]
                U_k = w_k_all[1]
                V_k = w_k_all[2]
                
                C_inf = w_inf_all[6]
                rho_inf = w_inf_all[0]
                P_inf = w_inf_all[4]
                U_inf = w_inf_all[1]
                V_inf = w_inf_all[2]
                
                vec_tau = np.array([dx, dy]) / ds
                vec_norm = np.array([dy, -dx]) / ds
                Un_k = np.dot(vec_norm, w_k_all[1:3])
                Ut_k = np.dot(vec_tau, w_k_all[1:3])
                Un_inf = np.dot(vec_norm, w_inf_all[1:3])
                Ut_inf = np.dot(vec_tau, w_inf_all[1:3])
                
                R_p = Un_k + 2 * C_k / (self.case.gamma - 1)
                R_m = Un_inf - 2 * C_inf / (self.case.gamma - 1)
                Un_edge = (R_p + R_m) / 2
                C_edge = (R_p - R_m) * (self.case.gamma - 1) / 4
                Man_inf = np.abs(Un_inf / C_inf)
                
                if Un_edge <= 0: # 入流 inflow
                    if Man_inf <= 1: # 亚声速 入流 subsonic inflow
                        s = P_inf / rho_inf ** self.case.gamma
                        rho_edge = (C_edge ** 2 / s / self.case.gamma) ** (1 / (self.case.gamma-1))
                        P_edge = s * rho_edge ** self.case.gamma
                        U_edge = U_inf + (Un_edge - Un_inf) * vec_tau[1]
                        V_edge = V_inf + (Un_edge - Un_inf) * vec_tau[0]
                        Z = Un_edge * ds
                        flux = np.array([
                            rho_edge * Z,
                            Z*rho_edge*U_edge + P_edge * dy,
                            Z*rho_edge*V_edge - P_edge * dx,
                            Z*rho_edge*( self.case.gamma * P_edge / (self.case.gamma - 1) / rho_edge + (U_edge**2 + V_edge**2) / 2)
                        ])
                        Q[k] = Q[k] + flux
                    else:
                        # 超声速 入流 supersonic inflow
                        rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.case.w_inf_all
                        flux = self.case.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                        Q[k] = Q[k] + flux
                        pass
                else: # 出流
                    if Man_inf <= 1: # 亚声速 出流 subsonic outflow
                        s = P_k / rho_k ** self.case.gamma
                        rho_edge = (C_edge ** 2 / s / self.case.gamma) ** (1 / (self.case.gamma-1))
                        P_edge = s * rho_edge ** self.case.gamma
                        U_edge = U_k + (Un_edge - Un_k) * vec_tau[1]
                        V_edge = V_k + (Un_edge - Un_k) * vec_tau[0]
                        Z = Un_edge * ds
                        flux = np.array([
                            rho_edge * Z,
                            Z*rho_edge*U_edge + P_edge * dy,
                            Z*rho_edge*V_edge - P_edge * dx,
                            Z*rho_edge*( self.case.gamma * P_edge / (self.case.gamma - 1) / rho_edge + (U_edge**2 + V_edge**2) / 2)
                        ])
                        Q[k] = Q[k] + flux
                        pass
                    else:
                        # 超声速 出流 supersonic outflow
                        rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.W_cur_all[k]
                        flux = self.case.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                        Q[k] = Q[k] + flux
                        pass
                    pass
                pass
            else:
                # 通量计算
                w_k = self.W_cur[k]
                w_p = self.W_cur[p]
                w_edge = (w_k + w_p) / 2
                rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.case.conserve_to_all(w_edge)
                flux = self.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                Q[k] = Q[k] + flux
                Q[p] = Q[p] - flux
                # 粘性计算
                # 二阶粘性
                alpha_edge = np.abs( U_edge * dy - V_edge * dx ) + C_edge * ds
                P_p = self.W_cur_all[p][4]
                P_k = self.W_cur_all[k][4]
                nu_edge = np.abs( (P_p - P_k) / (P_p + P_k) )
                epsilon_edge2 = nu_edge * self.case.k2
                d2_edge = alpha_edge * epsilon_edge2 * (w_p - w_k)
                D2[k] = D2[k] + d2_edge
                D2[p] = D2[p] - d2_edge
                # 四阶粘性
                epsilon_edge4 = np.max([0, self.case.k4 - epsilon_edge2])
                d4_edge = -alpha_edge * epsilon_edge4 * (self.W_cur_laplace[p] - self.W_cur_laplace[k])
                D4[k] = D4[k] + d4_edge
                D4[p] = D4[p] - d4_edge
                
                t[k] = t[k] + alpha_edge
                t[p] = t[p] + alpha_edge
                pass
            pass
        t = self.case.CFL * self.grid.vol / t
        self.Q = Q
        self.D = D2 + D4
        self.t = t
        pass
    
    def each_step_non_vis(self):
        self.conserve_to_all()
        self.calculate_laplace()
        Q = np.zeros([self.grid.ncells, 4])
        for edge in range(self.grid.nedges):
            i, j, k, p = self.grid.iedge[edge]
            dx = self.grid.dx[edge]
            dy = self.grid.dy[edge]
            ds = self.grid.ds[edge]
            if p == -2:
                w_k = self.W_cur[k]
                p = self.W_cur_all[k][4]
                flux = np.array([
                    0, p * dy, -p * dx, 0
                ])
                Q[k] = Q[k] + flux
                pass
            elif p == -3:
                w_k = self.W_cur[k]
                w_k_all = self.W_cur_all[k]
                w_inf = self.case.w_inf
                w_inf_all = self.case.w_inf_all
                
                C_k = w_k_all[6]
                rho_k = w_k_all[0]
                P_k = w_k_all[4]
                U_k = w_k_all[1]
                V_k = w_k_all[2]
                
                C_inf = w_inf_all[6]
                rho_inf = w_inf_all[0]
                P_inf = w_inf_all[4]
                U_inf = w_inf_all[1]
                V_inf = w_inf_all[2]
                
                vec_tau = np.array([dx, dy]) / ds
                vec_norm = np.array([dy, -dx]) / ds
                Un_k = np.dot(vec_norm, w_k_all[1:3])
                Ut_k = np.dot(vec_tau, w_k_all[1:3])
                Un_inf = np.dot(vec_norm, w_inf_all[1:3])
                Ut_inf = np.dot(vec_tau, w_inf_all[1:3])
                
                R_p = Un_k + 2 * C_k / (self.case.gamma - 1)
                R_m = Un_inf - 2 * C_inf / (self.case.gamma - 1)
                Un_edge = (R_p + R_m) / 2
                C_edge = (R_p - R_m) * (self.case.gamma - 1) / 4
                Man_inf = np.abs(Un_inf / C_inf)
                
                if Un_edge <= 0: # 入流 inflow
                    if Man_inf <= 1: # 亚声速 入流 subsonic inflow
                        s = P_inf / rho_inf ** self.case.gamma
                        rho_edge = (C_edge ** 2 / s / self.case.gamma) ** (1 / (self.case.gamma-1))
                        P_edge = s * rho_edge ** self.case.gamma
                        U_edge = U_inf + (Un_edge - Un_inf) * vec_tau[1]
                        V_edge = V_inf + (Un_edge - Un_inf) * vec_tau[0]
                        Z = Un_edge * ds
                        flux = np.array([
                            rho_edge * Z,
                            Z*rho_edge*U_edge + P_edge * dy,
                            Z*rho_edge*V_edge - P_edge * dx,
                            Z*rho_edge*( self.case.gamma * P_edge / (self.case.gamma - 1) / rho_edge + (U_edge**2 + V_edge**2) / 2)
                        ])
                        Q[k] = Q[k] + flux
                    else:
                        # 超声速 入流 supersonic inflow
                        rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.case.w_inf_all
                        flux = self.case.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                        Q[k] = Q[k] + flux
                        pass
                else: # 出流
                    if Man_inf <= 1: # 亚声速 出流 subsonic outflow
                        s = P_k / rho_k ** self.case.gamma
                        rho_edge = (C_edge ** 2 / s / self.case.gamma) ** (1 / (self.case.gamma-1))
                        P_edge = s * rho_edge ** self.case.gamma
                        U_edge = U_k + (Un_edge - Un_k) * vec_tau[1]
                        V_edge = V_k + (Un_edge - Un_k) * vec_tau[0]
                        Z = Un_edge * ds
                        flux = np.array([
                            rho_edge * Z,
                            Z*rho_edge*U_edge + P_edge * dy,
                            Z*rho_edge*V_edge - P_edge * dx,
                            Z*rho_edge*( self.case.gamma * P_edge / (self.case.gamma - 1) / rho_edge + (U_edge**2 + V_edge**2) / 2)
                        ])
                        Q[k] = Q[k] + flux
                        pass
                    else:
                        # 超声速 出流 supersonic outflow
                        rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.W_cur_all[k]
                        flux = self.case.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                        Q[k] = Q[k] + flux
                        pass
                    pass
                pass
            else:
                # 通量计算
                w_k = self.W_cur[k]
                w_p = self.W_cur[p]
                w_edge = (w_k + w_p) / 2
                rho_edge, U_edge, V_edge, E_edge, P_edge, H_edge, C_edge = self.case.conserve_to_all(w_edge)
                flux = self.calculate_flux(rho_edge, U_edge, V_edge, P_edge, H_edge, dx, dy)
                Q[k] = Q[k] + flux
                Q[p] = Q[p] - flux
                pass
            pass
        self.Q = Q
        pass
    
    def rk4_iteration(self):
        self.each_step() # 需要计算粘性通量的
        dt = np.min(self.t)
        W0 = copy.copy(self.W_cur)
        for step in range(4):
            self.W_cur = W0 - ((dt * rk4_coeff[step] * (self.Q - self.D)).T / self.grid.vol).T
            self.each_step_non_vis()
            pass
        residual = np.max(np.abs(self.W_cur - W0))
        return residual / self.residual_base
        pass
    
    def post_process(self):
        rho_node = np.zeros(self.grid.nnodes)
        vol_add_node = np.zeros(self.grid.nnodes)
        for k in range(self.grid.ncells):
            for i in range(3):
                ji = self.grid.icell[k, i]
                rho_node[ji] = rho_node[ji] + self.W_cur[k][0] * self.grid.vol[k]
                vol_add_node[ji] = vol_add_node[ji] + self.grid.vol[k]
                pass
            pass
        self.rho_result = rho_node / vol_add_node
        pass
    
    def start_simulate(self):
        self.initialize_field()
        self.residual_base = 1
        self.rk4_iteration()
        self.residual = np.zeros(self.case.STEP)
        self.residual[0] = 1
        self.residual_base = np.max(np.abs(self.W_cur - self.case.w_inf))
        for step in range(1, self.case.STEP):
            self.residual[step] = self.rk4_iteration()
            print("-------------------------------------")
            print("step: " + str(step))
            print("residual: " +  str(self.residual[step]))
            print("-------------------------------------")
            pass
        self.post_process()
        pass
    
    pass

###################################################################################################