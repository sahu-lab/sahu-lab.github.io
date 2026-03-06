### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 53a7e44d-c430-4eba-b7a6-43c593ed2dea
md"""
# Q2. Ideal gas work in the Carnot cycle $$\textcolor{red}{\texttt{[14]}}$$

In this problem, we are concerned with using numerical integration to calculate the work done in a [Carnot cycle](https://en.wikipedia.org/wiki/Carnot_cycle), where the working fluid is an ideal gas.
You will learn about the importance of the Carnot cycle in your thermodynamics course.
However, we are concerned with using numerical integration to tell us about the work done in a cycle.

To begin, suppose our system consists of 1 mole of helium gas, which can be approximated as ideal.
You are all familiar with the ideal gas equation of state, which---in the context of this problem---can be understood to represent the pressure as a function of the volume:
```math
P
\, = \, \dfrac{N R \, T}{V}
~.
```
The Carnot cycle is illustrated below, and consists of four steps between the four states (the latter of which are labeled 1, 2, 3, and 4).

![](https://upload.wikimedia.org/wikipedia/commons/0/06/Carnot_cycle_p-V_diagram.svg)

- __Isothermal expansion__ at a "hot" temperature `T_hot` (labeled as $$T_1$$ in the figure).
  Since the temperature is constant, we have
  ```math
  P_{1 \rightarrow 2}
  \, = \, \dfrac{N R \, T_{\text{hot}}}{V}
  ~.
  ```

- __Isentropic expansion,__ where there is no heat transfer.
  As you will learn from the study of thermodynamics, the relation between pressure and volume can be expressed as
  ```math
  P_{2 \rightarrow 3} \cdot \big( V_{2 \rightarrow 3} \big)^\gamma
  \, = \, \text{constant}
  ~,
  ```
  where $$\gamma$$ is the [heat capacity ratio](https://en.wikipedia.org/wiki/Heat_capacity_ratio)---which, for simple monatomic gases such as helium, is known to be
  ```math
  \gamma
  \, = \, \dfrac{5}{3}
  ~.
  ```

- __Isothermal compression__ at a "cold" temperature `T_cold` (labeled as $$T_2$$ in the figure), for which
  ```math
  P_{3 \rightarrow 4}
  \, = \, \dfrac{N R \, T_{\text{cold}}}{V}
  ~.
  ```

- __Isentropic compression,__ where there is no heat transfer.
  Once again, we have
  ```math
  P_{4 \rightarrow 1} \cdot \big( V_{4 \rightarrow 1} \big)^\gamma
  \, = \, \text{constant}
  ~.
  ```

The total work of the Carnot cycle can be calculated as
```math
\int_{\text{cycle}} P \, \text{d} V
~.
```
You will use numerical integration to calculate this work.
"""

# ╔═╡ f55127ef-6635-4057-9e83-7ad877a68bce
md"""
## Numerical integration $$\textcolor{red}{\texttt{[4]}}$$

Write a function
> `calc_integral(f::Function, a::Number, b::Number, N::Number)`
that numerically calculates
```math
\int_a^b f(x) \, \text{d} x
```
by splitting up the integral into `N` elements, and then approximating the area of each element.
You may use any of the integration methods we discussed in lecture.
Make sure that you correctly handle scenarios where $$b < a$$, for which
```math
\int_a^b f(x) \, \text{d} x
\, = \, - \int_b^a f(x) \, \text{d} x
~.
```

"""

# ╔═╡ 8ba043dc-ee13-4c4a-960f-9fb30b05d753
md"""
## Work calculation $$\textcolor{red}{\texttt{[10]}}$$

Write a function
> `carnot_work(N::Number)`
that calculates the total work done over one cycle, in Joules.
Note that this work is the area of the four-sided region surrounded by blue curves in the above figure.
Here, `N` denotes the number of elements that each integral is split up into.
The following constants are defined below:
- `T_cold` and `T_hot` are the cold and hot temperatures, in Kelvin
- `P_3` is the pressure at point 3 in the above figure, in Pascals
- `P_4` is the pressure at point 4 in the above figure, in Pascals (it must be larger than `P_3`)
- `gamma` is dimensionless, and its value is provided above
- `R` is the ideal gas constant, expressed in Joules / (mole * Kelvin)

To calculate the total work, the function `carnot_work` should call `calc_integral` four times.
Much of the work in writing this function is ensuring the proper limits of integration.

> __*Feel free to define any additional functions, as necessary.
> Make sure your code will still work if we change any of the below quantities (do not hard code anything to these specific values)*__


"""

# ╔═╡ 24b61c6c-28c2-46ff-8b06-c613c687d329
T_hot = 600.0

# ╔═╡ 97d22bda-fe51-44ef-999b-6b284034162b
T_cold = 300.0

# ╔═╡ 5222623d-e357-4ba7-ab52-d0256a16f010
P_3 = 1e5

# ╔═╡ 08c88e28-1769-4b51-b0ed-fb8705f3bd94
P_4 = 2e5

# ╔═╡ db88c52c-ba60-46da-9a0d-9e05d383d32c
gamma = 5/3

# ╔═╡ 703076a6-6939-4af3-9d5d-2199eb824d53
R = 8.314

# ╔═╡ 167b5fe8-3bdd-430c-95d5-616f287b39fc
md"""
> Evaluate your function for `N = 10, 33, 100, 333, 1000`
"""

# ╔═╡ Cell order:
# ╟─53a7e44d-c430-4eba-b7a6-43c593ed2dea
# ╟─f55127ef-6635-4057-9e83-7ad877a68bce
# ╟─8ba043dc-ee13-4c4a-960f-9fb30b05d753
# ╠═24b61c6c-28c2-46ff-8b06-c613c687d329
# ╠═97d22bda-fe51-44ef-999b-6b284034162b
# ╠═5222623d-e357-4ba7-ab52-d0256a16f010
# ╠═08c88e28-1769-4b51-b0ed-fb8705f3bd94
# ╠═db88c52c-ba60-46da-9a0d-9e05d383d32c
# ╠═703076a6-6939-4af3-9d5d-2199eb824d53
# ╟─167b5fe8-3bdd-430c-95d5-616f287b39fc
