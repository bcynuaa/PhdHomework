import sympy as spy
import numpy as np
import matplotlib.pyplot as plt

L = spy.symbols("L")
E = spy.symbols("E")
q = spy.symbols("q")
J1 = spy.symbols("J_1")
J2 = spy.symbols("J_2")
ML = spy.symbols("M_L")
x = spy.symbols("x")
EJ1 = E*J1
EJ2 = E*J2

theta_M = ML*L / EJ1 + ML*L / EJ2
print(spy.latex(theta_M))
theta_M

theta_q1 = -q*L*L**2 / 2 / EJ1 - (q*L**2/2)*L / EJ1 
print(spy.latex(theta_q1))
theta_q1

theta_q2 = -q*L**3 / 6 / EJ2
print(spy.latex(theta_q2))
theta_q2

theta_q = theta_q1 + theta_q2
print(spy.latex(theta_q))
theta_q

eq = spy.Eq(theta_M+theta_q, 0)
print(spy.latex(eq))
eq

ML_res = spy.solve(eq, ML)[0]
print(spy.latex(ML_res))
ML_res

theta_M_x_1 = ML_res*x / EJ1
print(spy.latex(theta_M_x_1))
theta_M_x_1

theta_q_x_1 = -q*L * x**2 / 2 / EJ1 - (q*L*(L-x+L/2))*x / EJ1
theta_q_x_1 = theta_q_x_1.simplify()
print(spy.latex(theta_q_x_1))
theta_q_x_1

theta_x_1 = theta_M_x_1 + theta_q_x_1
print(spy.latex(theta_x_1))
theta_x_1

theta_end_1 = theta_x_1.subs(x, L)
print(spy.latex(theta_end_1))
theta_end_1

theta_M_x_2 = ML_res*(x-L) / EJ2
print(spy.latex(theta_M_x_2))
theta_M_x_2

theta_q_x_2 = -q*(x-L)**3 / 6 / EJ2 - q*(2*L-x)*(x-L)**2 / 2 / EJ2 - (q*(2*L-x)**2/2)*(x-L) / EJ2
theta_q_x_2 = theta_q_x_2.simplify()
print(spy.latex(theta_q_x_2))
theta_q_x_2

theta_x_2 = theta_M_x_2 + theta_q_x_2 + theta_end_1
theta_x_2 = theta_x_2.simplify()
print(spy.latex(theta_x_2))
theta_x_2

theta_x_2.subs(x, 2*L).simplify()

theta_x = spy.Piecewise((theta_x_1, x <= L), (theta_x_2, x > L))
theta_x = theta_x.simplify()
print(spy.latex(theta_x))
theta_x

w_M_x_1 = ML_res*x**2 / 2 / EJ1
print(spy.latex(w_M_x_1))
w_M_x_1

w_q_x_1 = -q*L * x**3 / 3 / EJ1 - (q*L*(L-x+L/2))*x**2 / 2 / EJ1
w_q_x_1 = w_q_x_1.simplify()
print(spy.latex(w_q_x_1))
w_q_x_1

w_x_1 = w_M_x_1 + w_q_x_1
w_x_1 = w_x_1.simplify()
print(spy.latex(w_x_1))
w_x_1

(w_x_1.diff(x) - theta_x_1).simplify()

w_end_1 = w_x_1.subs(x, L)
w_end_1 = w_end_1.simplify()
print(spy.latex(w_end_1))
w_end_1

w_M_x_2 = ML_res*(x-L)**2 / 2 / EJ2
w_M_x_2 = w_M_x_2.simplify()
print(spy.latex(w_M_x_2))
w_M_x_2

w_q_x_2 = -q*(x-L)**4 / 8 / EJ2 - q*(2*L-x)*(x-L)**3 / 3 / EJ2 - (q*(2*L-x)**2/2) * (x-L)**2 / 2 / EJ2
w_q_x_2 = w_q_x_2.simplify()
print(spy.latex(w_q_x_2))
w_q_x_2

w_x_2 = w_M_x_2 + w_q_x_2 + w_end_1 + theta_end_1*(x-L)
w_x_2 = w_x_2.expand()
print(spy.latex(w_x_2))
w_x_2

(w_x_2.diff(x).simplify() - theta_x_2).simplify()

w_x = spy.Piecewise((w_x_1, x <= L), (w_x_2, x > L))
w_x = w_x.simplify()
print(spy.latex(w_x))
w_x

subs_dict = {
    L: 0.15,
    E: 207e9,
    J1: (38.1e-3)**4*np.pi/64,
    J2: (50.8e-3)**4*np.pi/64,
    q: 35e3,
}

theta_x_numerical = theta_x.subs(subs_dict)
theta_x_numerical = theta_x_numerical.simplify()
print(spy.latex(theta_x_numerical))
theta_x_numerical

w_x_numerical = w_x.subs(subs_dict)
w_x_numerical = w_x_numerical.simplify()
print(spy.latex(w_x_numerical))
w_x_numerical

plt.figure(figsize=(10, 5), facecolor="white")
x_values = np.linspace(0, 0.3, 1000)
theta_x_values = np.array([theta_x_numerical.subs(x, x_value) for x_value in x_values])
plt.plot(x_values, theta_x_values)
plt.xlabel("$x$ [m]")
plt.ylabel(r"$\theta(x)$ [rad]")
plt.grid()
plt.title(r"Deflection angle $\theta(x)$")
plt.savefig("../images/deflection_angle_analytical.pdf", bbox_inches="tight")
plt.savefig("../images/deflection_angle_analytical.png", bbox_inches="tight")
plt.show()

plt.figure(figsize=(10, 5), facecolor="white")
x_values = np.linspace(0, 0.3, 1000)
w_x_values = np.array([w_x_numerical.subs(x, x_value) for x_value in x_values])
plt.plot(x_values, w_x_values)
plt.xlabel("$x$ [m]")
plt.ylabel(r"$w(x)$ [m]")
plt.grid()
plt.title(r"Deflection $w(x)$")
plt.savefig("../images/deflection_analytical.pdf", bbox_inches="tight")
plt.savefig("../images/deflection_analytical.png", bbox_inches="tight")
plt.show()

from sympy.utilities.codegen import codegen

# generate python code for theta(x)
theta_x_code = codegen(
    ("analyticalAngle", theta_x_numerical),
    "julia",
    "analyticalAngle",
    header=False,
    empty=False,
)
print(theta_x_code[0][1])

# generate python code for w(x)
w_x_code = codegen(
    ("analyticalDeflection", w_x_numerical),
    "julia",
    "analyticalDeflection",
    header=False,
    empty=False,
)
print(w_x_code[0][1])