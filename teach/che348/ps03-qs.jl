### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ fa166b10-5b3b-11ef-17ef-a3ec27649460
begin
	using LinearAlgebra
end

# ╔═╡ 689ca386-cac2-40d7-87e3-ec88ef576107
md"""
### Initializing packages
"""

# ╔═╡ 0dd37eac-dec0-43b8-a5ae-92e77328a572
md"""
# Q1: Eigenvalues of projections in $$\mathbb{R}^3$$ $$\textcolor{red}{\texttt{[3]}}$$

In Problem Set 2, you consider a specific example of a projection matrix $$\boldsymbol{P}$$ in $$\mathbb{R}^3$$ that would act on a vector and project it in the direction of another arbitrary vector $$\boldsymbol{u}$$.
You also verified the following relations, which generally hold:

-  $$\boldsymbol{P} \boldsymbol{u} - \boldsymbol{u} = \boldsymbol{0}$$


-  $$\boldsymbol{P}^2 - \boldsymbol{P} = \boldsymbol{0}$$


-  $$(\boldsymbol{a} - \boldsymbol{P} \boldsymbol{a}) \boldsymbol{\cdot} \boldsymbol{u} = 0 \qquad$$ for any vector $$\boldsymbol{a} \in \mathbb{R}^3$$


Without doing any calculations or writing any code, justify what the eigenvalues of $$\boldsymbol{P}$$ are.
$$\textcolor{red}{\texttt{[3]}}$$

"""

# ╔═╡ 1ddb5bdd-009d-4162-bd23-64d8f534d08a
md"""
[Type your answer here]
"""

# ╔═╡ 6a81ba5c-fa83-456c-bec7-fef8284a88a5
md"""
# Q2: Matrix inverse vs linear solve $$\textcolor{red}{\texttt{[2]}}$$

We return to the generic linear equation
```math
\boldsymbol{A}
\boldsymbol{y}
\, = \, \boldsymbol{b}
~,
```
where the matrix $$\boldsymbol{A}$$ and vector $$\boldsymbol{b}$$ are known.
We seek to determine the unknown vector $$\boldsymbol{y}$$, for which there are two options:

- first, determine $$\boldsymbol{A}^{-1}$$, and then multiply it with $$\boldsymbol{b}$$


- use the built in "backslash" command, "`\`"

To demonstrate that one of these is much faster than the other, we will randomly generate large vectors and matrices:

"""

# ╔═╡ 43175e03-e00c-4720-b2db-cdbc5ab018d5
size = 2000

# ╔═╡ 29734760-7b3d-42ca-8832-f85a614a3b78
A = rand(Float64, (size, size))

# ╔═╡ 32cd60c6-eebd-4417-93bb-20ab0ee51b89
b = rand(Float64, size)

# ╔═╡ af3ab563-f7f6-4e89-bd91-d40f1f1883a3
md"""
Now, use the `@time` command (as done in [chapter 2](https://sahu-lab.github.io/teaching/che348/ch02-linear-algebra.html)) to see whether `A^{-1} * b` or `A \ b` is faster.
Make sure to run the cell (by pressing `[Shift]` + `[Enter]`) several times, to get reliable data.
$$\textcolor{red}{\texttt{[2]}}$$

"""

# ╔═╡ edf58ea6-7287-4295-93fa-3c69d4421434
md"""
# Q3: Practice with loops & functions $$\textcolor{red}{\texttt{[2]}}$$

- write a function named `f1` that takes in a number `n` and returns `3*n + n^2`.
  Evaluate `f1(2)` and `f1(3)` to verify that your function is working $$\textcolor{red}{\texttt{[1]}}$$


- write a function named `f2` that takes in a number `n` (assumed to be a positive integer) and returns the sum of all even numbers less than or equal to `n`. Evaluate `f2(5)` and `f2(8)` to verify that your function is working $$\textcolor{red}{\texttt{[1]}}$$
"""

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
# ╟─689ca386-cac2-40d7-87e3-ec88ef576107
# ╠═fa166b10-5b3b-11ef-17ef-a3ec27649460
# ╟─0dd37eac-dec0-43b8-a5ae-92e77328a572
# ╠═1ddb5bdd-009d-4162-bd23-64d8f534d08a
# ╟─6a81ba5c-fa83-456c-bec7-fef8284a88a5
# ╠═43175e03-e00c-4720-b2db-cdbc5ab018d5
# ╠═29734760-7b3d-42ca-8832-f85a614a3b78
# ╠═32cd60c6-eebd-4417-93bb-20ab0ee51b89
# ╟─af3ab563-f7f6-4e89-bd91-d40f1f1883a3
# ╟─edf58ea6-7287-4295-93fa-3c69d4421434
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
