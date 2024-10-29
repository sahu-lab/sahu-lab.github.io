### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d675a8f6-9592-11ef-01b0-89d810b1964b
md"""
# Dimensional Analysis

Thus far, we have been remarkably lax about dimensions and units in our discussions.
Doing so allowed us to focus on numerical details as we got started.
However, in engineering practice (and in science more broadly) dimensions are one of the most important entities to keep track of!
We will therefore be mindful of units and dimensions from this point forward.
First, some terminology:
- __*dimensions*__ refer to fundamental physical characteristics, such as mass, length, time, and charge


- __*units*__ refer to human-made (but arbitrary) reference quantities, such as grams, meters, seconds, and Coulombs
You have most likely had some experience with converting between different units—for example, converting between meters and feet or kilograms and pounds.
To begin our more rigorous discussion, we present two principles that seem straightforward but contain deep truths.


"""

# ╔═╡ 1481bceb-717f-4a40-b86a-3372116c6b67
md"""
## Principle 1: Two quantities with different dimensions cannot be added or subtracted

- 2 seconds + 3 seconds = 5 seconds $\checkmark$


- 3 meters + 5 meters = 8 meters $\checkmark$


- 5 meters + 8 seconds = ???

> __*Note:*__ it is of course possible to add two quantities of different __units,__ as long as they have the same __dimension.__
> To do so, we require the appropriate *conversion factor* between units

!!! question "Example Question"
    What is 1 meter + 1 foot?

This principle is especially important when considering differential equations!
Consider, for example, the dynamics of a mass attached to a spring—which you studied in your physics courses.
The statement of Newton's second law ($$\boldsymbol{F} = m \boldsymbol{a}$$) in one dimension reads
```math
m \, \dfrac{\mathrm{d}^2 x}{\mathrm{d} t^2}
\, + \, k \, x
\, = \, 0
~.
```
Let's say you forgot the dimensions of the spring constant $$k$$.
You can determine it with dimensional analysis.
```math
m \, \dfrac{\mathrm{d}^2 x}{\mathrm{d} t^2}
\qquad
\text{ has dimensions of }
~
\text{mass} \ast \dfrac{\text{length}}{\text{time}^2}
```
```math
\Rightarrow \quad
k \, x
\qquad
\text{ must have dimensions of }
~
\text{mass} \ast \dfrac{\text{length}}{\text{time}^2}
```
```math
x
~\text{ has dimensions of length}
\qquad
\Rightarrow
\qquad
k
~ \text{ has dimensions of }
~ \dfrac{\text{mass}}{\text{time}^2}
```

With this information, let us try to determine (without solving the differential equation) the *characteristic time* of mass oscillations.
Since
```math
k
~ [=] \ \dfrac{\text{mass}}{\text{time}^2}
\qquad
\text{and}
\qquad
k
~ [=] \ \text{mass}
~,
```
we know
```math
\dfrac{m}{k}
~ [=] \ \text{time}^2
\qquad
\Rightarrow
\qquad
\sqrt{\dfrac{m}{k} \,}
~ [=] \ \text{time}
~.
```
Up to a factor of $$2 \pi$$, we found the period of oscillation without doing any calculations!


"""

# ╔═╡ 723dfadb-43a1-44ad-b7d7-1c9a2ac3dbe6
md"""
## Principle 2: Our choice of units cannot affect any physical outcome

- If you choose to express the gravitational constant as $$g = 9.8 ~ \text{m}/\text{sec}^2$$ and I express it as $$g = 32.2 ~ \text{ft}/\text{sec}^2$$, we __*must*__ get the same answer

!!! info "Important Idea"
	Nature does not care about what arbitrary units we choose

- By a similar logic, if we are referring to positions in space relative to some coordinate system, the results cannot depend on the origin of our coordinate system or the way we positions our axes!


!!! question "Conceptual Question"
	We have repeatedly considered the scalar differential equation
	```math
	\dfrac{\text{d} y}{\text{d} t}
	\, = \, k \, y
	\qquad
	\text{with}
	\quad
	y (t = 0)
	\, = \, y_0
	~.
	```
	What are the dimensions of $$k$$?
"""

# ╔═╡ b1d015b2-a908-490b-a527-695a85940a79
md"""
## A return to prior concepts

When first taking general chemistry, you may have been confused about the dimensions of different reaction-rate constants, depending on whether the reaction was first-, second-, or third-order.
It turns out that these units are simply required for dimensional consistency.
Consider the equation for the time evolution of the concentration of a single species, $A$, in a well-mixed batch reactor:
```math
\dfrac{\text{d} c_A}{\text{d} t}
\, = \, r^{}_A
~,
```
Where $r_A$ is the so-called reaction rate.
If $A$ is consumed in a first-order reaction, then
```math
\dfrac{\text{d} c_A}{\text{d} t}
\, = \, - k^{(1)} \, c_A
~.
```
If $A$ is produced in a second-order reaction involving another species $B$, then
```math
\dfrac{\text{d} c_A}{\text{d} t}
\, = \, k^{(2)} \, c_A \, c_B
~.
```
If $A$ is consumed in a third-order reaction, one possibility is for
```math
\dfrac{\text{d} c_A}{\text{d} t}
\, = \, k^{(3)} \, c^2_A \, c_B
~.
```
Assuming all concentrations are measured in moles per leter, and time is measured in seconds, what are the units of $k^{(1)}$, $k^{(2)}$, and $k^{(3)}$? 

"""

# ╔═╡ Cell order:
# ╟─d675a8f6-9592-11ef-01b0-89d810b1964b
# ╟─1481bceb-717f-4a40-b86a-3372116c6b67
# ╟─723dfadb-43a1-44ad-b7d7-1c9a2ac3dbe6
# ╟─b1d015b2-a908-490b-a527-695a85940a79
