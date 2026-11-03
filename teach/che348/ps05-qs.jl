### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 07d3e656-7437-11ef-3342-0bea283a9848
begin
	using LinearAlgebra
end

# ╔═╡ 53a7e44d-c430-4eba-b7a6-43c593ed2dea
md"""
# Q2. Gradients and Jacobians in $$\mathbb{R}^2$$ $$\textcolor{red}{\texttt{[14]}}$$

In this problem, we are concerned with quantities in $$\mathbb{R}^2$$.
However, you should think about how to extend these ideas and concepts to higher dimensions, with the computer taking care of repetetive steps.

Write a function `f(x)` that takes as input a $$2 \times 1$$ vector $$\boldsymbol{x}$$ and outputs the scalar $$\textcolor{red}{\texttt{[1]}}$$
```math
f(\boldsymbol{x})
\, = \, 4 - (x_1)^2 - 2 x_1 x_2 + 3 x_1 - (x_2)^2 + 2 x_2
~.
```

"""

# ╔═╡ 78458c59-d00b-4807-8459-a59c8739ba13
md"""
---
"""

# ╔═╡ b227266e-39ee-4953-91fc-bcb1380b8d7b
md"""
Write a function `df(x, h)` that numerically approximates $$\boldsymbol{\nabla} f$$, and returns it as a $$2 \times 1$$ vector.
$$\textcolor{red}{\texttt{[3]}}$$
"""

# ╔═╡ 4963a4db-4587-4587-aa47-45ac2ffbd9dd
md"""
---
"""

# ╔═╡ d143ae97-b96f-4e19-9318-84d746a8b0ff
md"""
Check that your function `df` calculates the correct gradient by evaluating it at the points $$(1, 2)$$ and $$(3, -1)$$.
Use $$h = 0.001$$.
In order to carry out this check, calculate the analytical gradient and write a function `df_exact(x)` which returns this quantity.
In each case, calculate the norm of the difference between the exact and approximate gradient.
$$\textcolor{red}{\texttt{[4]}}$$
"""

# ╔═╡ 2f3dc665-641c-4a8e-9827-edf099d6f8e0
md"""
---
"""

# ╔═╡ dfb733ce-b026-4eea-b157-407a6ff32ec8
md"""
Now consider the vector function $$\boldsymbol{g} (\boldsymbol{x})$$, where both $$\boldsymbol{g}$$ and $$\boldsymbol{x}$$ are elements of $$\mathbb{R}^2$$.
In particular,
```math
\boldsymbol{g} (\boldsymbol{x})
\, = \, \begin{bmatrix}
g_1 (x_1, x_2) \\
g_2 (x_1, x_2) \\
\end{bmatrix}
~,
```
where
```math
\begin{aligned}
g_1 (x_1, x_2)
\, &= \, 3 (x_1)^2 - 2 (x_2)^3 + 5
\\[5pt]
g_2 (x_1, x_2)
\, &= \, (x_1 + x_2)^2 - 3 x_1 (x_2)^2 + 2
~.
\end{aligned}
```
Write a function `g(x)` that accepts a $$2 \times 1$$ vector $$\boldsymbol{x}$$ as an argument, and outputs the $$2 \times 1$$ vector $$\boldsymbol{g} (\boldsymbol{x})$$. 
$$\textcolor{red}{\texttt{[2]}}$$
"""

# ╔═╡ 45fde833-594e-4ff9-98ae-70f234f8dd60
md"""
---
"""

# ╔═╡ fcd994ea-ad9d-4089-8709-2e8f10c3bfa0
md"""
Fill out the function `dg(x, h)` below, which accepts a $$2 \times 1$$ vector $$\boldsymbol{x}$$ as an argument and outputs the $$2 \times 2$$ Jacobian matrix
```math
\boldsymbol{\nabla} \boldsymbol{g}
\, = \, \begin{bmatrix}
\dfrac{\partial g_1}{\partial x_1} & \dfrac{\partial g_1}{\partial x_2} \\[5pt]
\dfrac{\partial g_2}{\partial x_1} & \dfrac{\partial g_2}{\partial x_2} \\[5pt]
\end{bmatrix}
```
The function `dg` should reference only the earlier function `g`, and calculate the Jacobian numerically.
$$\textcolor{red}{\texttt{[4]}}$$
"""

# ╔═╡ 943d8605-a48e-4ce1-bc1d-44515263c01c
function dg(x, h)
	return zeros(2, 2)
end

# ╔═╡ e6971f04-e702-4636-91ea-5d28cbff8428
md"""
If you implemented the function `dg(x, h)` correctly, then the code below should output something close to the analytical result
```math
\boldsymbol{\nabla} \boldsymbol{g} \big\rvert_{(1, \, 3)}
\, = \, \begin{bmatrix}
6 & -54 \\
-19 & -10 \\
\end{bmatrix}
~.
```
"""

# ╔═╡ 18048151-0003-4423-a1c3-0d263f0d0897
dg([1, 3], 0.001)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╠═07d3e656-7437-11ef-3342-0bea283a9848
# ╟─53a7e44d-c430-4eba-b7a6-43c593ed2dea
# ╟─78458c59-d00b-4807-8459-a59c8739ba13
# ╟─b227266e-39ee-4953-91fc-bcb1380b8d7b
# ╟─4963a4db-4587-4587-aa47-45ac2ffbd9dd
# ╟─d143ae97-b96f-4e19-9318-84d746a8b0ff
# ╟─2f3dc665-641c-4a8e-9827-edf099d6f8e0
# ╟─dfb733ce-b026-4eea-b157-407a6ff32ec8
# ╟─45fde833-594e-4ff9-98ae-70f234f8dd60
# ╟─fcd994ea-ad9d-4089-8709-2e8f10c3bfa0
# ╠═943d8605-a48e-4ce1-bc1d-44515263c01c
# ╟─e6971f04-e702-4636-91ea-5d28cbff8428
# ╠═18048151-0003-4423-a1c3-0d263f0d0897
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
