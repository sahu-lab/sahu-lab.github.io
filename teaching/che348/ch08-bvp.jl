### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ fa645e3a-a09c-11ef-11d1-7f5c746c9ada
md"""
# Boundary Value Problems

For the initial value problems studied previously, all conditions arose at time $t = 0$, and we used information at the current time to approximate the solution at some nearby time in the future (recall that $\Delta t$ must be small in order for our methods to apply).
We are now concerned with a different class of problems, which often arise in science and engineering, that involve spatial (rather than temporal) domains.
In such situations, one often applies boundary conditions at each end of the domain (in one dimension) or along the entire boundary (in multiple dimensions).
We will leave the study of equations involving multiple spatial dimensions when we learn how to solve partial differential equations.



### Example: Plug flow reactor

Let us consider a typical example: a [plug flow reactor](https://en.wikipedia.org/wiki/Plug_flow_reactor_model) that is operating at steady state.

![](https://upload.wikimedia.org/wikipedia/commons/3/3b/Pipe-PFR.svg)

Fluid containing a solute is flowed through a reactor of length $\ell$ and cross-sectional area $A$ with velocity $v$.
Due to a catalyist contained in the reactor, a chemical reaction proceeds with a known rate law $r$.
If there is a single chemical species with concentration $c$ and the rate law is first-order, then at steady state a mass balance yields
```math
v \, \dfrac{\text{d} c}{\text{d} z}
\, = \, D \, \dfrac{\text{d}^2 c}{\text{d} z^2}
\, - \, k c
~.
```
The above equation is fairly general, and is sometimes called the convection--diffusion--reaction equation, as it contains terms that describe each of these phenomena: the term on the left-hand side describes the convection of the chemical species due to the fluid velocity $v$, the first term on the right-hand side captures the diffusion of the chemical species with a diffusion constant $D$, and the last term tells us how the concentration changes due to the chemical reaction.
Since the equation contains two spatial derivatives, two boundary conditions are required, and in general one boundary condition is prescribed on each end of the domain.

!!! question "Important Question"
	What are the dimensions of all the quantities in the above equation?

In your kinetics course, you will likely study the scenario in which __convection dominates diffusion,__ for which
```math
v \, \dfrac{\text{d} c}{\text{d} z}
\, = \, - \, k c
~,
```
and only one boundary condition is required—for example, $c(0) = c_0$.

!!! question "Conceptual Question"
	How can the ideas we learned about __initial value problems__ be applied to this reduced scenario?


### Example: Catalyst pellet

As mentioned earlier, the convection--diffusion--reaction equation is fairly general, and arises in many situations.
Suppose we are interested in understanding the steady-state concentration profile inside a single catalyis pellet, which can be modeled as a thin, flat slab.
Within the pellet itself, there is no flow—and thus __no convection.__
The equation governing the 1-D concentration profile along the slab thickness is given by
```math
D \, \dfrac{\text{d}^2 c}{\text{d} z^2}
\, - \, k c
\, = \, 0
~,
```
which is sometimes called the reaction--diffusion equation.
You will see examples of these equations come up repeatedly in chemical engineering practice.


"""

# ╔═╡ 5ce4129a-d03b-448a-b7de-b662458147a0
md"""

## Finite difference method

While you are comfortable solving the above equations by hand, boundary-value problems can quickly become difficult to solve analytically.
We will require numerical tools to solve them.
Let us revisit the example of a catalyst pellet, which extends from $z = 0$ to $z = \ell$ (suppose it is much wider than it is tall).
Let there also be different concentrations on the two sides, for which
```math
D \, \dfrac{\text{d}^2 c}{\text{d} z^2}
\, - \, k c
\, = \, 0
\qquad
\text{for } ~
0 \le z \le \ell
~,
\qquad
\text{with } ~
c(0) = c^{}_0
~,
\quad
c(\ell) = c^{}_1
~.
```


### Non-dimensionalization

Our first step is to non-dimensionalize.
Following a different notation from before (you may use whichever you are more comfortable with), we define
```math
Z
\, \equiv \, \dfrac{z}{\ell}
\qquad
\text{and}
\qquad
C
\, \equiv \, \dfrac{c}{c^{}_0}
~.
```
We also define the dimensionless parameters
```math
\alpha
\, \equiv \, \dfrac{\ell^2 k}{D}
\qquad
\text{and}
\qquad
\beta
\, \equiv \, \dfrac{c^{}_1}{c^{}_0}
~,
```
and simplify our problem statement to
```math
\dfrac{\text{d}^2 C}{\text{d} Z^2}
\, - \, \alpha C
\, = \, 0
\qquad
\text{for } ~
0 \le Z \le 1
~,
\qquad
\text{with } ~
C(0) = 1
~,
\quad
C(1) = \beta
~.
```
At this point, there are only two parameters in the problem: $\alpha$, which compares the reaction and diffusion rates, and $\beta$, which captures the relative concentrations on the two sides of the slab.

!!! question "Practice Question"
	Determine the dimensionless form of the problem statement from the original one.


### Numerical solution

To solve for the concentration profile numerically, we first recognize there are an infinite number of positions, or $Z$-values, in the domain $Z \in [0, 1]$.
We cannot numerically satisfy the differential equation at an infinite number of points, and so we __discretize__ the domain into a finite number of intervals $N$.

|------------|------------|------------|------------|------------|------------|------------|------------|

If
```math
\Delta Z
\, := \, \dfrac{1}{N}
~,
```
is the width of a single interval, then the locations of the endpoints of each interval are $0$, $\Delta Z$, $2 \Delta Z$, $3\Delta Z$, $\ldots$, $N \Delta Z = 1$.
If $j$ is the interval number, then the endpoints of that interval are $(j - 1) \Delta Z$ and $ j \Delta Z$.
Our objective is to determine the concentrations at each of the interval endpoints, which are also called __nodes:__

> __Objective:__ $\text{find } ~ C(Z_j) \, , ~ \text{ where } ~ Z_j \equiv j \Delta Z ~ \text{ and } ~ j \in \{0, 1, 2, \ldots, N \}$

!!! info "Observation"
	There are $N + 1$ values of $j$, which means there are $N + 1$ unknowns.

Remember that two boundary conditions are provided: $C(0) = 1$ and $C(1) = \beta$.
In terms of our discretized (or approximate) concentration, this tells us
```math
C(Z_0)
\, = \, 1
\qquad
\text{and}
\qquad
C(Z_N)
\, = \, \beta
~.
```
Accordingly, we wish to determine the remaining $N-1$ unknowns: $C(Z_1)$, $C(Z_2)$, $\ldots$, $C(Z_{N-1})$.

!!! info "Shorthand"
	Let us introduce the notation
	```math
	C_j \, \equiv \, C(Z_j)
	```
	to reduce our writing.

For the $N-1$ unknowns, we require $N-1$ equations.
Recalling that the governing differential equation is valid everywhere on the domain, we choose to apply it at each of the $Z_j$ in the domain.
In other words,
```math
\dfrac{\text{d}^2 C}{\text{d} Z^2} \bigg\rvert_{Z_j}
\, - \, \alpha C_j
\, = \, 0
\qquad
\text{for } ~
j \in \{ 1, 2, \ldots, N-1 \}
~.
```
Our task now is to approximate the derivative numerically, given that the concentration $C$ is only known at a set of discrete points.
In an earlier problem set, you determined
```math
\dfrac{\text{d}^2 C}{\text{d} Z^2} \bigg\rvert_{Z_j}
\, \approx \, \dfrac{1}{(\Delta Z)^2} \, \Big(
C_{j+1}
\, - \, 2 C_j
\, + \, C_{j-1}
\Big)
~.
```
Upon substituting into the above expression, we find
```math
\dfrac{1}{(\Delta Z)^2} \, \Big(
C_{j+1}
\, - \, 2 C_j
\, + \, C_{j-1}
\Big)
\, - \, \alpha C_j
\, = \, 0
~.
```
With some rearrangement, we find
```math
C_{j+1}
\, - \, \big[
2 + \alpha (\Delta Z)^2
\big] C_j
\, + \, C_{j-1}
\, = \, 0
~.
```
To reduce our writing even further, we define the shorthand
```math
\gamma
\, \equiv \, 2 + \alpha (\Delta Z)^2
```
For example, when $j=1$, we have
```math
C_2
\, - \, \gamma \, C_1
\, + \, C_0
\, = \, 0
~,
```
where $C_0 = 1$ from the boundary condition.
When $j = 2$, we have
```math
C_3
\, - \, \gamma \, C_2
\, + \, C_1
\, = \, 0
~.
```
Notice that both equations involve $C_1$ and $C_2$.
In fact, $C_2$ also appears in the equation for $j = 3$:
```math
C_4
\, - \, \gamma \, C_3
\, + \, C_2
\, = \, 0
~.
```
Each equation, corresponding to a specific $j$, involves the concentrations at $j - 1$ and $j + 1$ as well.
Thus we cannot individually solve for each $C_j$ one at a time.
Instead, we need to solve for all of them together.




"""

# ╔═╡ Cell order:
# ╟─fa645e3a-a09c-11ef-11d1-7f5c746c9ada
# ╟─5ce4129a-d03b-448a-b7de-b662458147a0
