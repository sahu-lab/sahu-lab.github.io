### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 53a7e44d-c430-4eba-b7a6-43c593ed2dea
md"""
# Q2. Van der Waals fluid (continued) $$\textcolor{red}{\texttt{[18]}}$$

In the written portion of this problem set, you investigated some properties of the van der Waals fluid.
We will continue to analyze some of its properties, but now with numerical methods.
In all cases, we consider the dimensionless form of the equation, where now the '$$*$$' accent is dropped for notational convenience:
```math
p
\, = \, \dfrac{T}{v - 1}
\, - \, \dfrac{1}{v^2}
~.
```

"""

# ╔═╡ adafdfc7-d884-4860-a024-eebdcad30e96
md"""
## Pressure as a function $$\textcolor{red}{\texttt{[1]}}$$

Write a function `p_vdw(vol::Number, temp::Number)` that returns the (dimensionless) pressure of the van der Waals fluid, for a given (dimensionless) volume and (dimensionless) temperature.
For the remainder of this question, it is understood that all quantities are dimensionless.

"""

# ╔═╡ d851949b-f3c2-49af-a402-c5406ccbe56e
md"""
## Plotting the pressure $$\textcolor{red}{\texttt{[2]}}$$

Create a single plot of the pressure as a function of volume, which contains three curves (corresponding to three different temperatures).
One temperature should be the critical temperature $$T_{\text{c}}$$, one should be slightly above $$T_{\text{c}}$$, and one should be slighlty below $$T_{\text{c}}$$.
Remember, you calculated the critical temperature in the written portion of the problem set.


"""

# ╔═╡ 5b63121e-ea5b-4f1d-a6bf-9185e05839ab
md"""
## Numerical differentiation $$\textcolor{red}{\texttt{[2]}}$$

Write a function that numerically calculates the derivative $$\partial p / \partial v$$ of the van der Waals fluid.
You will need to choose
- the name of the function
- the arguments the function takes
- the numerical parameter $$h$$
Note that you should be able to call your function at different volumes and temperatures.

> Do **_NOT_** calculate $$\partial p / \partial v$$ analytically (by hand) here!

"""

# ╔═╡ f4c82f98-fe84-47a7-92f1-6d57879ce4e5
md"""
## Plotting the derivative $$\textcolor{red}{\texttt{[2]}}$$

At the three temperatures you chose above, plot $$\partial p / \partial v$$ using the function you wrote above.
All the curves should be in a single plot.
What qualitative difference do you see in the derivative, between temperatures above $$T_{\text{c}}$$ and temperatures below $$T_{\text{c}}$$?

"""

# ╔═╡ 23212f79-0fca-44ee-baba-587f8a733137
md"""
## Finding where the derivative is zero $$\textcolor{red}{\texttt{[5]}}$$

Adapt the bisection method we wrote in class to numerically determine where $$\partial p / \partial v = 0$$, for a given (i.e. user-specified) temperature.
Recall that for temperatures below $$T_{\text{c}}$$, there are two such volumes.
Your function should return both, in an array.
These are called **_spinodal points._**
"""

# ╔═╡ c284bd64-c687-40a0-81d0-9b7443052813
md"""
## Plotting the spinodal points $$\textcolor{red}{\texttt{[6]}}$$

In the exercise above, you considered a single temperature and found two volumes corresponding to spinodal points.
Now, sweep through many temperatures ranging between $$T = 0.11$$ and $$T = T_{\text{c}}$$ to generate the so-called __spinodal curve__ in the pressure--volume plane.
In plotting the spinodal curve, you do not need to note the temperature at which the values were found.
An example can be seen in the black, dash-dotted line in the figure below.

> Note that the pressure may in fact go negative, which is okay because the spinodal regions are unphysical.
> Additionally, this plot uses a different non-dimensionalization, so your results should not quantitatively agree with it.

![](https://upload.wikimedia.org/wikipedia/commons/6/6e/VdW_isotherms%2B2log.png)
"""

# ╔═╡ Cell order:
# ╟─53a7e44d-c430-4eba-b7a6-43c593ed2dea
# ╟─adafdfc7-d884-4860-a024-eebdcad30e96
# ╟─d851949b-f3c2-49af-a402-c5406ccbe56e
# ╟─5b63121e-ea5b-4f1d-a6bf-9185e05839ab
# ╟─f4c82f98-fe84-47a7-92f1-6d57879ce4e5
# ╟─23212f79-0fca-44ee-baba-587f8a733137
# ╟─c284bd64-c687-40a0-81d0-9b7443052813
