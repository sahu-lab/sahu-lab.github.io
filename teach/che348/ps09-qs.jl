### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5811cf9b-25ef-4851-9d51-8c7ec233cd03
begin
	using PlutoUI
	using LinearAlgebra
end

# ╔═╡ d433673a-9af0-11ef-31f5-532e0f539822
md"""
# Q2. Dynamics of a chemical reactor $$\textcolor{red}{\texttt{[14]}}$$

Here, we will continue analyzing the scenario described in the written portion of the problem set.
Recall that in terms of the dimensionless concentration and temperature, we have the coupled equations
```math
\dfrac{\text{d} c^*}{\text{d} t^*}
\, = \, 1
\, - \, c^*
\, - \, \tilde{K} \, (c^*)^2 \, e^{-\tilde{G}/T^*}
```
and
```math
\dfrac{\text{d} T^*}{\text{d} t^*}
\, = \, 1
\, + \, \tilde{U} \, T^*_\infty
\, - \, (1 + \tilde{U}) \, T^*
\, + \, \tilde{K} \, \tilde{H} \, (c^*)^2 \, e^{-\tilde{G}/T^*}
~.
```
In order to be able to solve these equations numerically, we want to express them in a form convenient to a computer.
To this end, we define
```math
y_1
\, \equiv \, c^*
\qquad
\text{and}
\qquad
y_2
\, \equiv \, T^*
~,
```
such that our collective unknowns are contained in the vector
```math
\boldsymbol{y}(t)
\, = \, \begin{bmatrix}
y_1 (t)
\\[3pt]
y_2 (t)
\end{bmatrix}
~.
```
Note that in our code, the "time" variable is understood to be the dimensionless time $t^*$.


"""

# ╔═╡ fcad97ac-2227-417b-9633-205f30b115bf
md"""
First, let us define the relevant constants in the problem.
We will drop the "$\sim$" accent for notational simplicity.
You can use the sliders to change the values of the different constants.
The slider value is printed for your convenience.

"""

# ╔═╡ 4389a709-c98b-4012-9642-c8f5fc2b24f4
@bind K Slider(1:10, default=3)

# ╔═╡ d22c5b6a-66e7-49a1-8efa-47995a81df0b
K

# ╔═╡ 8de396f8-1883-45fb-af70-33eb913fb5e1
@bind G Slider(0.1:0.1:2.0, default=0.1)

# ╔═╡ 0a3ae873-4b8b-4255-ae8c-bd220a9a4a42
G

# ╔═╡ 8022011a-7c65-479c-a879-f3d5472299f0
@bind U Slider(0.1:0.1:2.0, default=0.3)

# ╔═╡ de630fac-0faa-472e-b62a-cae1516a7812
U

# ╔═╡ ff2a0cdf-0718-4518-8e8d-4173337acd6f
@bind H Slider(-1.0:0.1:1.0, default=0.2)

# ╔═╡ 2ca31dab-a7fb-4003-ab1b-64d821f7ed89
H

# ╔═╡ d8f55bd9-7438-4b23-a332-2b2474e15ac1
@bind T_inf Slider(0.8:0.1:1.2, default=0.9)

# ╔═╡ fe4172e4-7e92-4315-8ef5-7e492741e3c6
T_inf

# ╔═╡ 613b487b-e836-44fa-80eb-72c3191b4c3c
md"""

## RHS function $$\textcolor{red}{\texttt{[3]}}$$

Write a function `rhs(y::Vector{Number})` that takes in the vector $\boldsymbol{y}$ and returns the right-hand side of the two differential equations provided above, collected into a vector.
Your function should use the values from the sliders above, i.e. `K`, `G`, `H`, `U`, and `T_inf`.

"""

# ╔═╡ 01bdab11-249e-48e5-b49d-482d2bb36c0a
function rhs(y)
	# the following return statement is incorrect
	return [y[1], y[2]]
end

# ╔═╡ 91c7dac4-e6cb-4cb8-92a7-3a61a67da47d
rhs([0.5, 0.8])[2]

# ╔═╡ 5e7c0923-8c43-4488-b55f-35449c135e92
md"""
!!! info "Check Your Work!"
	If your `rhs` function is implemented properly, then using the default values provided above you should find
	`rhs([0.5, 0.8]) = [-0.161873, 0.362375]`, as checked below
"""

# ╔═╡ 78196f94-7015-4b68-8cba-4d457753f192
if norm(rhs([0.5, 0.8]) - [-0.16187267693844665, 0.3623745353876893]) < 1e-10
	println("Your function is implemented correctly!")
else
	println("It looks like there's a mistake, keep at it!")
end

# ╔═╡ d22d34e2-afb8-4dff-b788-af0a73ed5a64
md"""
We are now interested in solving for the steady-states of the highly nonlinear functions above, using numerical methods.
In particular, we would like to find the vector `y_ss` such that `rhs(y_ss) = 0`.

"""

# ╔═╡ 0573a9fd-511d-4803-a831-f4e7774ee9d5
md"""
## Determining the steady-states numerically

Using the functions defined below, which correctly implement fixed-point iteration and the Newton--Raphson method, solve for `y_ss`.


"""

# ╔═╡ 16b2e39e-b2af-4eb4-82fe-0ac571765650
md"""

### Fixed-point iteration $$\textcolor{red}{\texttt{[5]}}$$

The function `fixed_pt(f, y0)` defined below solves for the vector `y` where `f(y) = y`.
However, we are interested in determining where `rhs(y) = 0`.
Thus, we need to define a __*new function*__ that can be passed into the `fixed_pt` function.
Recalling your written problem set from last week, write a function `fp_rhs(y)` for which `y = fp_rhs(y)` is a solution to `rhs(y) = 0`.
There are many options to choose from, but not all of them will lead to a converged solution.

"""

# ╔═╡ 4f292f8e-2746-41fb-b73a-49daa8b7d3d9
# Use fixed point iteration to find where f(y) = y, given an initial guess y0
function fixed_pt(f::Function, y0)
	eps_tol   = 1e-11
	max_count = 400

	y = y0
	count = 0
	error = norm(f(y) - y)
	while error > eps_tol && count < max_count
		y = f(y)
		error = norm(f(y) - y)
		count = count + 1
		# to be used for debugging
		#println(y, "\t", error)
	end

	# if the iteration reached the maximum count, or the values approach infinity
	# (thus leading to the error being "Not a Number", or NaN, then we failed)
	if count == max_count || isequal(error, NaN)
		println("fixed point not found!")
		return 0.0 * y0;
	else
		println("converged in ", count, " iterations!")
		return y
	end
end

# ╔═╡ 2416ac47-b542-453e-8486-230deeb508e2
md"""
Determine the steady-state value, according to the fixed-point method, and save it in a variable titled `y_ss_fp`.
"""

# ╔═╡ 82af52d8-a85b-481f-b282-3d071f84d961
md"""
### Newton--Raphson method $$\textcolor{red}{\texttt{[2]}}$$

Now determine the same steady-state solution using the function `newton_raphson` defined below.
This function takes a function `f` and an initial guess `y0`, and solves for `y` such that `f(y) = 0`.

"""

# ╔═╡ 3831d0b7-d8b8-4bc9-a578-8de2a3133428
# this function was adapted from the coding solution to problem set 5
function jacobian(f::Function, y)
	h = 1e-7
	num_dimensions = length(y)
	J = zeros(num_dimensions, num_dimensions)
	for i = 1:num_dimensions
		dy = zeros(num_dimensions)
		dy[i] = h
		J[:,i] = (f(y + dy) - f(y)) / h
	end

	return J
end

# ╔═╡ d254d9e6-c3cd-4ddf-94a4-50150353d4d8
# Use Newton--Raphson iteration to find where f(y) = 0, given an initial guess y0
function newton_raphson(f::Function, y0)
	eps_tol   = 1e-11
	max_count = 10

	y = y0
	count = 0
	error = norm(f(y))
	
	while error > eps_tol && count < max_count
		y = y - jacobian(f, y) \ f(y)
		error = norm(f(y))
		count = count + 1
		# to be used for debugging
		#println(y, "\t", error)
	end

	# if the iteration reached the maximum count, or the values approach infinity
	# (thus leading to the error being "Not a Number", or NaN, then we failed)
	if count == max_count || isequal(error, NaN)
		println("fixed point not found!")
		return 0.0 * y0;
	else
		println("converged in ", count, " iterations!")
		return y
	end
	
end

# ╔═╡ ce22e51c-85e3-4fe3-acb4-89fa3a6ef129
md"""
## Numerical stability of the steady state $$\textcolor{red}{\texttt{[4]}}$$

Now that we have determined the steady state solution, let us investigate the stability numerically.
Call the `jacobian` function defined above, using either your fixed-point or Newton--Raphson solution.

"""

# ╔═╡ 52b3d09c-0403-4d74-a4ed-0efd296b4447
md"""
!!! question "Question"
	Is the steady state stable or unstable?
	Why?
"""

# ╔═╡ e5893786-2a3b-4ac0-b918-d67d90d7766d
md"""
[Type your answer here]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "33626088eef1949506fed17cbfc9d01571b54bb6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═5811cf9b-25ef-4851-9d51-8c7ec233cd03
# ╟─d433673a-9af0-11ef-31f5-532e0f539822
# ╟─fcad97ac-2227-417b-9633-205f30b115bf
# ╠═4389a709-c98b-4012-9642-c8f5fc2b24f4
# ╠═d22c5b6a-66e7-49a1-8efa-47995a81df0b
# ╠═8de396f8-1883-45fb-af70-33eb913fb5e1
# ╠═0a3ae873-4b8b-4255-ae8c-bd220a9a4a42
# ╠═8022011a-7c65-479c-a879-f3d5472299f0
# ╠═de630fac-0faa-472e-b62a-cae1516a7812
# ╠═ff2a0cdf-0718-4518-8e8d-4173337acd6f
# ╠═2ca31dab-a7fb-4003-ab1b-64d821f7ed89
# ╠═d8f55bd9-7438-4b23-a332-2b2474e15ac1
# ╠═fe4172e4-7e92-4315-8ef5-7e492741e3c6
# ╟─613b487b-e836-44fa-80eb-72c3191b4c3c
# ╠═01bdab11-249e-48e5-b49d-482d2bb36c0a
# ╠═91c7dac4-e6cb-4cb8-92a7-3a61a67da47d
# ╟─5e7c0923-8c43-4488-b55f-35449c135e92
# ╠═78196f94-7015-4b68-8cba-4d457753f192
# ╟─d22d34e2-afb8-4dff-b788-af0a73ed5a64
# ╟─0573a9fd-511d-4803-a831-f4e7774ee9d5
# ╟─16b2e39e-b2af-4eb4-82fe-0ac571765650
# ╠═4f292f8e-2746-41fb-b73a-49daa8b7d3d9
# ╟─2416ac47-b542-453e-8486-230deeb508e2
# ╟─82af52d8-a85b-481f-b282-3d071f84d961
# ╠═3831d0b7-d8b8-4bc9-a578-8de2a3133428
# ╠═d254d9e6-c3cd-4ddf-94a4-50150353d4d8
# ╟─ce22e51c-85e3-4fe3-acb4-89fa3a6ef129
# ╟─52b3d09c-0403-4d74-a4ed-0efd296b4447
# ╠═e5893786-2a3b-4ac0-b918-d67d90d7766d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
