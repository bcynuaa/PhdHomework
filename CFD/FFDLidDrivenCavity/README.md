# FFD_LidDrivenCavity

A repository containing FFD code to solve lid-driven cavity using FFD method.

The code was written in December 2022 when I studying FFD. Now I try to make it more readable and understandable.

What's more, the code repository plays the role as appendix for the 'Linear Algebra in System Control' course's final report.

# Requirements

- Language: Julia
- Pacakges:
    - `LinearAlgebra`
    - `SparseArrays`
    - `Plots`
    - `PyPlot`
    - `DelimitedFiles`
    - `CSV`
    - `DataFrames`

# Directory

- `src`: source code
- `example`: example of the code for $Re=100,400$ and grid number for $16\times 16, 32\times 32, 64\times 64$
- `image`: image of the contour and stream line, including vivid gif and the drawing `jupyter-notebook`
- `data`: data of the velocity at the middle of the cavity

# Brief Results

## $Re=100$

<figure class="half">
<center>
<img src="./image/Re100/u.gif" width=200/><img src="./image/Re100/v.gif" width=200/>

gif of the contour for $u$ and $v$
</center>
</figure

<figure class="half">
<center>
<img src="./image/Re100/streamline.png" width=200/><img src="./image/Re100/u_middle.png" width=200/>

stream line and $u$ at the middle of the cavity
</center>
</figure>

## $Re=400$

<figure class="half">
<center>
<img src="./image/Re400/u.gif" width=200/><img src="./image/Re400/v.gif" width=200/>

gif of the contour for $u$ and $v$
</center>
</figure>

<figure class="half">
<center>
<img src="./image/Re400/streamline.png" width=200/><img src="./image/Re400/u_middle.png" width=200/>

stream line and $u$ at the middle of the cavity
</center>
</figure>

# Reference

- [Fast Fluid Dynamics-people.sc.fsu.edu](https://people.sc.fsu.edu/~lb13f/projects/finite_difference/fast_fluid_dynamics.php#:~:text=Fast%20Fluid%20Dynamics%20%28FFD%29%20is%20a%20technique%20for,extended%20to%20wind%20load%20optimization%20among%20other%20applications.)
- [Github-Balavarun5/fast_fluid_dynamics](https://github.com/Balavarun5/fast_fluid_dynamics)
- [Nvidia-FFD](https://developer.nvidia.com/gpugems/gpugems/part-vi-beyond-triangles/chapter-38-fast-fluid-dynamics-simulation-gpu)
- [University of Colorado Boulder-FFD](https://www.colorado.edu/lab/sbs/fastfluiddynamics)