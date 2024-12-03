### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ f1ddb783-8fa2-424e-8a86-37577f4ddeef
md"""
# Dimensional Analysis

Thus far, we have been remarkably lax about dimensions and units in the problems we have solved.
Doing so allowed us to focus on the details of the numerics, and their derivations.
We have now reached the point where we would like to solve problems with real-world applications.
In order to do so, we must be mindful of dimensions and units.
As an added bonus, dimensional analysis often allows us to do less work when interpreting and analyzing scenarios numerically.


"""

# ╔═╡ 697bd288-e86c-45b1-8338-21b56d0eb326
md"""
## Chemical engineering example

Let us begin with a motivational example: the start-up of a [continuous stirred-tank reactor(CSTR)](https://en.wikipedia.org/wiki/Continuous_stirred-tank_reactor):

![](https://upload.wikimedia.org/wikipedia/commons/e/e0/Cstr.png)

For times $t < 0$, the reactor (with volume $V$) has an initial concentration of species $B$, $c^0_B$, and none of species $A$.
Then, at time $t = 0$, an inlet stream with flow rate $Q$ and species $A$ at concentration $c^{\text{in}}_{\! A}$ is started; an outleat stream exits with the same flow rate, and concentrations are at their well-mixed values.
In the reactor, suppose the reaction $A \rightarrow B$ occurs as a second-order reaction: $r^{}_{\! A} = - r^{}_{\! B} = - k \, c^2_{\! A}$.
The reactor is then governed by the system of equations
```math
\begin{aligned}
V \, \dfrac{\text{d} c^{}_{\! A}}{\text{d} t}
\, &= \, Q \, \big( c^{\text{in}}_{\! A} \, - \, c^{}_{\! A} \big)
\, - \, V \, k \, c^2_{\! A}
\\[8pt]
V \, \dfrac{\text{d} c^{}_{\! B}}{\text{d} t}
\, &= \, - Q \, c^{}_{\! B}
\, + \, V \, k \, c^2_{\! A}
~,
\end{aligned}
```
along with the initial conditions
```math
c^{}_{\! A} (0)
\, = \, 0
\qquad
\text{and}
\qquad
c^{}_{\! B} (0)
\, = \, c^0_{\! B}
~.
```
Thus, in modeling this scenario, there are two unknowns we need to solve for: $c^{}_{\! A} (t)$ and $c^{}_{\! B} (t)$.
In addition, the solution depends on five parameters: $V$, $Q$, $c^{\text{in}}_{\! A}$, $c^0_{\! B}$, and $k$.

!!! question "Important Question"
	Do all five parameters uniquely affect the final solution?

	__*ABSOLUTELY NOT!*__

We recognize that to solve the problem, we should set *characteristic scales* for quantities with different dimensions.
To this end, we define
```math
t^*
\, = \, \dfrac{t}{\tau}
~,
\qquad
c^*_{\! A}
\, = \, \dfrac{c^{}_{\! A}}{c^{}_{\text{ref}}}
~,
\qquad
\text{and}
\qquad
c^*_{\! B}
\, = \, \dfrac{c^{}_{\! B}}{c^{}_{\text{ref}}}
~,
```
where thus far we have not specified either the characteristic timescale $\tau$ or the reference concentration $c^{}_{\text{ref}}$.
By substituting the dimensionless quantities into the governing equation for $c^{}_{\! A}$, we find
```math
\dfrac{V c^{}_{\text{ref}}}{\tau} \, \dfrac{\text{d} c^{*}_{\! A}}{\text{d} t^*}
\, = \, Q \, \big( c^{\text{in}}_{\! A} \, - \, c^{}_{\text{ref}} \, c^{*}_{\! A} \big)
\, - \, V \, k \, c^{2}_{\text{ref}} \, ( c^*_{\! A} )^2
~.
```
After some rearrangement, we obtain
```math
\dfrac{\text{d} c^{*}_{\! A}}{\text{d} t^*}
\, = \, \dfrac{\tau \, Q}{V} \, \bigg( \dfrac{c^{\text{in}}_{\! A}}{c^{}_{\text{ref}}} \, - \, c^{*}_{\! A} \bigg)
\, - \, \tau \, k \, c^{}_{\text{ref}} \, ( c^*_{\! A} )^2
~.
```
At this point, we choose for
```math
\tau
\, = \, \dfrac{V}{Q}
\qquad
\text{and}
\qquad
c^{}_{\text{ref}}
\, = \, c^{\text{in}}_{\! A}
~,
```
and additionally define the dimensionless rate constant
```math
k^*
\, = \, \dfrac{V \, k \, c^{\text{in}}_{\! A}}{Q}
~.
```
With this, the dynamical equation for the concentration of species $A$ is presented in dimensionless form as
```math
\dfrac{\text{d} c^{*}_{\! A}}{\text{d} t^*}
\, = \, 1 \, - \, c^{*}_{\! A}
\, - \, k^* ( c^*_{\! A} )^2
~.
```
After some algebra, we similarly find the equation for the concentration of species $B$ to be given by
```math
\dfrac{\text{d} c^{*}_{\! B}}{\text{d} t^*}
\, = \, - \, c^{*}_{\! B}
\, + \, k^* ( c^*_{\! A} )^2
~.
```
Finally, the initial conditions are expressed as
```math
c^*_{\! A} (0)
\, = \, 0
\qquad
\text{and}
\qquad
c^*_{\! B} (0)
\, = \, \dfrac{c^{0}_{\! B}}{c^{\text{in}}_{\! A}}
```
We thus find there are only two independent parameters: $k^*$ and $c^{0}_{\! B} / c^{\text{in}}_{\! A}$.
Dimensional analysis has greatly simplified the problem!

Our task now is to recast the problem into one that we can solve numerically.
To this end, we define
```math
y_1
\, \equiv \, c^*_{\! A}
~,
\qquad
y_2
\, \equiv \, c^*_{\! B}
~,
\qquad
\text{and}
\qquad
\boldsymbol{y}
\, = \, \begin{bmatrix}
y_1
\\[4pt]
y_2
\end{bmatrix}
~.
```
In addition, we define
```math
f_1 (y_1, y_2)
\, \equiv \, 1
\, - \, y_1
\, - \, k^* (y_1)^2
~,
\qquad
f_2 (y_1, y_2)
\, \equiv \, - \, y_2
\, + \, k^* (y_1)^2
~,
```
and
```math
\boldsymbol{f} (\boldsymbol{y})
\, = \, \begin{bmatrix}
f_1(y_1, y_2)
\\[4pt]
f_2(y_1, y_2)
\end{bmatrix}
~.
```
We thus express the system of equations in the general form
```math
\dfrac{\text{d} \boldsymbol{y}}{\text{d} t}
\, = \, \boldsymbol{f} (\boldsymbol{y})
~.
```


"""

# ╔═╡ d675a8f6-9592-11ef-01b0-89d810b1964b
md"""
## Basic ideas

Now that we have seen the utility of dimensional analysis, we present the most important ideas.
First, some terminology:
- __*dimensions*__ refer to fundamental physical characteristics, such as mass, length, time, and charge


- __*units*__ refer to human-made (but arbitrary) reference quantities, such as grams, meters, seconds, and Coulombs
You most likely have experience converting between different units—for example, converting between meters and feet or kilograms and pounds.
We begin a more rigorous discussion by presenting two important principles.


"""

# ╔═╡ 1481bceb-717f-4a40-b86a-3372116c6b67
md"""
### Principle 1: Dimensional equality

Two quantities can be added or subtracted only if they have the same dimensions.
Thus, in any equation—including in differential equations—all terms must have the same dimensions.
For the simplest example:

- 2 seconds + 3 seconds = 5 seconds $\checkmark$


- 3 meters + 5 meters = 8 meters $\checkmark$


- 5 meters + 8 seconds = ???

> __*Note:*__ it is of course possible to add two quantities of different __units,__ as long as they have the same __dimension.__
> To do so, we require the appropriate *conversion factor* between units


### Principle 2: Invariance to units

While dimensions are fundamental to nature, units are human-made.
Accordingly, our choice of units cannot affect any of our physical predictions.
For example:

- If you choose to express the gravitational constant as $$g = 9.8 ~ \text{m}/\text{sec}^2$$ and I express it as $$g = 32.2 ~ \text{ft}/\text{sec}^2$$, we __*must*__ get the same answer


- By a similar logic, if we are referring to positions in space relative to some coordinate system, the results cannot depend on the origin of our coordinate system or the way we positions our axes


## Additional examples

The subject of dimensional analysis is complex, and we don't have time to go through concepts in great detail.
Instead, we will look at a couple more examples that will set us up for our solution of differential equations.


### Mass on a spring

Consider the dynamics of a mass attached to a spring—which you studied in your physics courses.
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
~~ [=] ~~ \text{mass} \ast \dfrac{\text{length}}{\text{time}^2}
```
```math
\Rightarrow \quad
k \, x
~~ [=] ~~ \text{mass} \ast \dfrac{\text{length}}{\text{time}^2}
```
```math
x
~~ [=] ~~ ~\text{length}
\qquad
\Rightarrow
\qquad
k
~~ [=] ~~ \dfrac{\text{mass}}{\text{time}^2}
~.
```

With this information, let us try to determine (without solving the differential equation) the *characteristic time* of mass oscillations.
Since
```math
k
~~ [=] ~~ \dfrac{\text{mass}}{\text{time}^2}
\qquad
\text{and}
\qquad
k
~~ [=] ~~ \text{mass}
~,
```
we know
```math
\dfrac{m}{k}
~~ [=] ~~ \text{time}^2
\qquad
\Rightarrow
\qquad
\sqrt{\dfrac{m}{k} \,}
~~ [=] ~~ \text{time}
~.
```
Up to a factor of $$2 \pi$$, we found the period of oscillation without doing any calculations!


"""

# ╔═╡ b1d015b2-a908-490b-a527-695a85940a79
md"""
### Dimensions of reaction-rate constants

When first taking general chemistry, you may have been confused about the dimensions of different reaction-rate constants, depending on whether the reaction was first-, second-, or third-order.
It turns out that these units are simply required for dimensional consistency.
Consider the equation for the time evolution of the concentration of a single species, $A$, in a well-mixed batch reactor:
```math
\dfrac{\text{d} c^{}_{\! A}}{\text{d} t}
\, = \, r^{}_{\! A}
~,
```
Where $r^{}_{\! A}$ is the so-called reaction rate.
If $A$ is consumed in a first-order reaction, then
```math
\dfrac{\text{d} c^{}_{\! A}}{\text{d} t}
\, = \, - k^{(1)} \, c^{}_{\! A}
~.
```
If $A$ is produced in a second-order reaction involving another species $B$, then
```math
\dfrac{\text{d} c^{}_{\! A}}{\text{d} t}
\, = \, k^{(2)} \, c^{}_{\! A} \, c^{}_{\! B}
~.
```
If $A$ is consumed in a third-order reaction, one possibility is for
```math
\dfrac{\text{d} c^{}_{\! A}}{\text{d} t}
\, = \, k^{(3)} \, c^{2}_{\! A} \, c^{}_{\! B}
~.
```
Assuming all concentrations are measured in moles per leter, and time is measured in seconds, what are the units of $k^{(1)}$, $k^{(2)}$, and $k^{(3)}$? 

"""

# ╔═╡ 3d9b55e6-9868-4300-bad6-3110021caa41
md"""
### Spring–mass–dampener

Now let us suppose that in the spring–mass system, there is some friction that acts against the motion of the mass.
The simplest mathematical description of such a force is one proportional to the velocity, for which
```math
m \, \dfrac{\mathrm{d}^2 x}{\mathrm{d} t^2}
\, + \, \mu \, \dfrac{\mathrm{d} x}{\mathrm{d} t}
\, + \, k \, x
\, = \, 0
~.
```

!!! question "Dimensional analysis question"
	What are the dimensions of the drag coefficient $\mu$?

For a second-order differential equation, we require two boundary conditions.
Let us choose
```math
x(0)
\, = \, x_0
\qquad
\text{and}
\qquad
\dfrac{\mathrm{d} x}{\mathrm{d} t} (0)
\, = \, 0
~,
```
where $x_0$ is a known constant.
We are interested in solving for $x(t)$, given four parameters: $m$, $\mu$, $k$, and $x_0$.
We can once again ask if these four parameters independently affect the solution.

Let us non-dimensionalize the problem by defining
```math
x^*
\, = \, \dfrac{x}{x_0}
\qquad
\text{and}
\qquad
t^*
\, = \, \dfrac{t}{\tau}
~,
```
where for now $\tau$ is an unknown time scale.
The boundary conditions simplify to
```math
x^* (0)
\, = \, 1
\qquad
\text{and}
\qquad
\dfrac{\mathrm{d} x^*}{\mathrm{d} t^*} (0)
\, = \, 0
~.
```
Moreover, the governing differential equation can be expressed as
```math
\dfrac{m \, x_0}{\tau^2} \, \dfrac{\mathrm{d}^2 x^*}{\mathrm{d} (t^*)^2}
\, + \, \dfrac{\mu \, x_0}{\tau} \, \dfrac{\mathrm{d} x^*}{\mathrm{d} t^*}
\, + \, k \, x_0 \, x^*
\, = \, 0
~.
```
With some rearrangement, we find
```math
\dfrac{\mathrm{d}^2 x^*}{\mathrm{d} (t^*)^2}
\, + \, \dfrac{\mu \, \tau}{m} \, \dfrac{\mathrm{d} x^*}{\mathrm{d} t^*}
\, + \, \dfrac{k \, \tau^2}{m} \, x^*
\, = \, 0
~.
```
From here, there are two natural choices for the yet-to-be-specified time scale $\tau$:
```math
\tau^{}_1
\, = \, \dfrac{m}{\mu}
\qquad
\text{and}
\qquad
\tau^{}_2
\, = \, \sqrt{\dfrac{m}{k} \,}
~.
```
The first timescale is the decay time if there was no spring, while the second timescale is the one we found earlier when considering the period of oscillations in the absence of any drag.
We choose
```math
\tau
\, = \, \tau^{}_2
\, = \, \sqrt{\dfrac{m}{k} \,}
~,
```
for which our differential equation simplifies to
```math
\dfrac{\mathrm{d}^2 x^*}{\mathrm{d} (t^*)^2}
\, + \, \dfrac{\mu}{\sqrt{k \, m \,}} \, \dfrac{\mathrm{d} x^*}{\mathrm{d} t^*}
\, + \, x^*
\, = \, 0
~.
```
Finally, for notational simplicity, we define
```math
\alpha
\, \equiv \, \dfrac{\mu}{\sqrt{k \, m \,}}
```
as the ratio of timescales.
Our governing equation and boundary conditions are then written as
```math
\dfrac{\mathrm{d}^2 x^*}{\mathrm{d} (t^*)^2}
\, + \, \alpha \, \dfrac{\mathrm{d} x^*}{\mathrm{d} t^*}
\, + \, x^*
\, = \, 0
~,
```
with
```math
x^* (0)
\, = \, 1
\qquad
\text{and}
\qquad
\dfrac{\mathrm{d} x^*}{\mathrm{d} t^*} (0)
\, = \, 0
~.
```
We have thus found there is only a single relevant dimensionless parameter in this problem, as opposed the original four.

As we will see in the next chapter, our numerical implementaitons will generally require the problem statement to be expressed as $\boldsymbol{\dot{y}} = \boldsymbol{f} (\boldsymbol{y})$.
To this end, we define
```math
y_1
\, \equiv \, x^*
~,
\qquad
y_2
\, \equiv \, \dfrac{\mathrm{d} x^*}{\mathrm{d} t^*}
~,
\qquad
f_1 (y_1, y_2)
\, \equiv \, y_2
~,
\qquad
f_2 (y_1, y_2)
\, \equiv \, - y_1
\, - \, \alpha y_2
~,
```
along with
```math
\boldsymbol{y}
\, = \, \begin{bmatrix}
y_1
\\[4pt]
y_2
\end{bmatrix}
\qquad
\text{and}
\qquad
\boldsymbol{f} (\boldsymbol{y})
\, = \, \begin{bmatrix}
f_1 (y_1, y_2)
\\[4pt]
f_2 (y_1, y_2)
\end{bmatrix}
~.
```
The dynamics of the spring--mass--dampener system are then indeed given by
```math
\boldsymbol{\dot{y}}
\, = \, \boldsymbol{f} (\boldsymbol{y})
~.
```

"""

# ╔═╡ Cell order:
# ╟─f1ddb783-8fa2-424e-8a86-37577f4ddeef
# ╟─697bd288-e86c-45b1-8338-21b56d0eb326
# ╟─d675a8f6-9592-11ef-01b0-89d810b1964b
# ╟─1481bceb-717f-4a40-b86a-3372116c6b67
# ╟─b1d015b2-a908-490b-a527-695a85940a79
# ╟─3d9b55e6-9868-4300-bad6-3110021caa41
