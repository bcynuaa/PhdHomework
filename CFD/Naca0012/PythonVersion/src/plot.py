# import os
# os.chdir("..//src")
from project import *

grid_file = "..//data//naca0012.grd"
grid = Grid(grid_file)
case = Case(0.8, 2.5)
solver = Solver(case, grid)
solver.start_simulate()

triangle = tri.Triangulation(solver.grid.xy.T[0], solver.grid.xy.T[1], solver.grid.icell)

plt.figure(figsize=(6, 4), dpi=200, facecolor="white")
plt.gca().set_aspect(1)

plt.tricontourf(triangle, solver.rho_result, cmap="jet")
plt.colorbar()

plt.xlim(-1.5, 2.5)
plt.ylim(-2, 2)

plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title(r"$\rho$")
plt.savefig("..//image//rho_contour.png", bbox_inches="tight")
plt.savefig("..//image//rho_contour.pdf", bbox_inches="tight")
plt.show()