### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 1f42dd9e-a82f-498d-a6b2-73d1a05603fd
begin
	using Plots
end

# ╔═╡ 8268c060-58b4-11ef-0d04-dd88ba82f850
md"""

# Algebraic equations

Earlier in the course, we considered equations of the form
```math
\boldsymbol{A} \, \boldsymbol{y}
\, = \, \boldsymbol{b}
~,
```
where
$$\boldsymbol{y} \in \mathbb{R}^n$$
and
$$\boldsymbol{b} \in \mathbb{R}^n$$,
are $$n$$-dimensional vectors, and
$$\boldsymbol{A} \in \mathbb{R}^{n \times n}$$
is a square $$n \times n$$ matrix.
__*Linear equations*__ of this form are (in principle) straightforward to solve: assuming $$\boldsymbol{A}$$ is invertible, we determine the unknown vector $$\boldsymbol{y}$$ as
```math
\boldsymbol{y}
\, = \, \boldsymbol{A}^{-1} \boldsymbol{b}
~.
```

We now focus on the solution of __*nonlinear equations,*__ which is significantly more complicated.
Notationally, we seek a solution $$\boldsymbol{y}$$ of the equation
```math
\boldsymbol{f} (\boldsymbol{y})
\, = \, \boldsymbol{0}
~,
```
where
$$\boldsymbol{f} \in \mathbb{R}^n$$
and
$$\boldsymbol{y} \in \mathbb{R}^n$$.
The basic concepts will be demonstrated for the simplest case, where
$$n = 1$$
and we look to solve equations of the form
```math
f(y)
\, = \, 0
~.
```
To this end, consider the following:

!!! question "Motivating Question"
	Using arithmetic operators alone ($$+$$, $$-$$, $$\times$$, $$\div$$), how would you determine the numerical value of $$\sqrt{3}$$ to within a given tolerance $$\epsilon_{\text{tol}}$$?
	In other words, how would you determine a value $$v$$ such that
	```math
	\big\lvert 3 - v^2 \big\rvert \, < \, \epsilon_{\text{tol}}
	~,
	```
	where (for example)
	$$\epsilon_{\text{tol}} = 10^{-12}$$.


"""

# ╔═╡ a349ab6c-6d43-4850-9848-21d1d2158882
sqrt(3)

# ╔═╡ 1ffa3210-1442-499c-a066-1a9bb174f6a4
√3

# ╔═╡ e37983bb-04e5-4d01-bb15-9700774014c1
3 - (√3)^2

# ╔═╡ aca8795a-a924-4718-b266-aec2a648a049
md"""

## Bisection method

For the specific case of determining the numerical value of $$\sqrt{3}$$, we know that
$$\sqrt{1} = 1$$
and
$$\sqrt{4} = 2$$,
and so
```math
1 \, < \, \sqrt{3} \, < \, 2
~.
```
To put this on a firm mathematical footing, we introduce the function
```math
f(y)
\, = \, y^2
\, - \, 3
~,
```
and seek to determine $$\bar{y}$$ where
$$f(\bar{y}) = 0$$.

"""

# ╔═╡ 9b937cc8-f70b-476f-b194-18bd5ca7fc91
f(y) = y^2 - 3

# ╔═╡ 151957bb-cb42-43e8-9d9f-240b1776eaaf
f(1)

# ╔═╡ 2df83b27-7dcf-4293-b390-c8434bb2a8f9
f(2)

# ╔═╡ ac0791ee-059b-4f03-b984-78bc90a2667a
md"""

Since
$$f(1) = -2 < 0$$
and
$$f(2) = 1 > 0$$,
the [intermediate value theorem](https://en.wikipedia.org/wiki/Intermediate_value_theorem) guarantees
```math
\bar{y}
\, \in \, (1, 2)
~.
```
Let us thus choose our __*initial guess*__ $$y_0$$ for $$\bar{y}$$ as the midpoint:
$$y_0 = 1.5$$,
which __*bisects*__ the domain.

"""

# ╔═╡ 3028cddf-d42f-49af-9753-57f0dd1b55eb
y0 = (1 + 2)/2

# ╔═╡ fe19ad8f-6275-4248-9e70-e3a8d6cf0bde
md"""

Next, we evaluate the function $$f$$ at the initial guess $$y_0$$:

"""

# ╔═╡ 9fb354dd-e0b3-4a49-b81d-af91e679c799
f(y0)

# ╔═╡ de00f330-a035-4adb-bc88-3cd1c55edf75
md"""

Since
$$f(y_0) < 0$$
and
$$f(2) > 0$$,
we now know that
$$\bar{y} \in (y_0, 2)$$.
Thus, as our next guess, we choose $$y_1$$ to be the midpoint of the new, smaller domain:

"""

# ╔═╡ 44c90555-ce2d-435b-93ad-ae58bee105c4
y1 = (y0 + 2)/2

# ╔═╡ 72a4f34e-47dd-4be3-aaf7-ace5bd47f166
md"""

At this point, we repeat the process:

"""

# ╔═╡ 24ef0a11-86cd-44d9-8695-4fa5480b9b0e
f(y1) # greater than 0

# ╔═╡ 5dc003d7-bea8-4744-9f0f-6a16151dcd9e
y2 = (y0 + y1) / 2 # bisect the new domain, from y0 to y1

# ╔═╡ f474de65-1416-4e7f-9855-5bc920e7314a
f(y2) # find f(y2) < 0, so the new domain is (y2, y1)

# ╔═╡ f036f742-e5c7-45d2-91a6-3ccb0fe890be
md"""

In this manner, we generate a sequence of approximations to $$\bar{y}$$:
```math
\big\{ \, y_0, \, y_1, \, y_2, \, \ldots \, \big\}
~,
\qquad
\text{where}
\qquad
\lim_{n \rightarrow \infty} y_n
\, = \, \bar{y}
~.
```
Since each iteration bisects the length of the interval in which we are looking, this method is called the [bisection method](https://en.wikipedia.org/wiki/Bisection_method).

!!! question "Practical Question"
	Recall $$\sqrt{3}$$ is __*irrational,*__ and so it contains infinitely many digits.
	It is thus __*impossible*__ for us to represent its __*exact*__ value in a computer, as there is a maximum number of digits we can use!
	> *How can we calculate a number $$y_n$$ that is within the specified tolerance, i.e. for which*
	> 
	> $$\big\lvert \sqrt{3} \, - \, y_n \big\rvert \, < \, \epsilon_{\text{tol}} ~?$$
	
	Also, remember the concept of __*machine precision:*__ there is a limit to the accuracy of *any calculation we carry out!* 
"""

# ╔═╡ 98540be7-1c08-4d39-9ecc-1d9bac36c913
eps()

# ╔═╡ f7571334-7267-46e3-9418-1fb5aea32f54
md"""

### Error analysis

In the example above, we chose the inital domain for $$y$$ to be from $$1 \rightarrow 2$$.
Let us denote the length of this domain, in the general case, as $$\ell$$.
Then, since $$y_0$$ is the midpoint of this range, the __*maximum error*__ between our choice of $$y_0$$ and the true solution $$\bar{y}$$ is $$\ell / 2$$.
If we call the error at any iteration $$n$$ as $$r_n$$, then
```math
r_0
\, \le \, \dfrac{\ell}{2}
\qquad
\text{where}
\qquad
r_n
\, := \, \big\lvert \, y_n \, - \, \bar{y} \, \big\rvert
\quad
\text{for}
\quad
n = 0, 1, 2, \ldots
~.
```
The domains on either side of $$y_0$$ have length $$\ell / 2$$.
When we choose $$y_1$$ to bisect one of these domains, we find
```math
r_1
\, \le \, \dfrac{\ell}{4}
\, = \, \dfrac{\ell}{2^2}
~.
```
By repeating this procedure, we find
```math
r_n
\, \le \, \dfrac{\ell}{2^{n+1}}
~.
```
Thus, to guarantee $$r_n \le \epsilon_{\text{tol}}$$, we solve for the smallest value of $$n$$ satisfying
```math
\dfrac{\ell}{2^{n+1}} \, \le \, \epsilon_{\text{tol}}
\qquad
\Rightarrow
\qquad
\dfrac{\ell}{\epsilon_{\text{tol}}} \, \le \, 2^{n+1}
\qquad
\Rightarrow
\qquad
n \, \ge \, \log_2 \bigg(\dfrac{\ell}{\epsilon_{\text{tol}}}\bigg) \, - \, 1
~.
```
Moreover, since $$n$$ is a natural number (it is the number of bisections), we can write
```math
n \, = \, \bigg\lceil
\log_2 \bigg(\dfrac{\ell}{\epsilon_{\text{tol}}}\bigg)
\bigg\rceil
\, - \, 1
~,
```
where we use the [ceiling function](https://en.wikipedia.org/wiki/Floor_and_ceiling_functions).



### Implementation

> We will code the numerical implementation together in class, as it is the most complex code we have written so far!



### Pros and cons

- upper bound on the error ✔
- intuitive reasoning ✔
- always works for continuous functions ✔
 
- does not work when $$\boldsymbol{y}$$ is a vector! ✖

"""

# ╔═╡ f2f1a24d-7e4c-4bda-ade3-101c455d2e51
md"""

## Fixed-point iteration

Recall that our objective is to find a solution of the equation
```math
0
\, = \, y^2
\, - \, 3
\, =: \, f(y)
~.
```
Any solution $$\bar{y}$$ to the above equation will also satisfy
$$0 = 3 - y^2$$,
and thus
```math
y
\, = \, 3
\, + \, y
\, - \, y^2
~.
```
Let us define
```math
g(y)
\, \equiv \, 3
\, + \, y
\, - \, y^2
~,
```
for which we seek a solution to equation
```math
y
\, = \, g(y)
~.
```
In general, such a solution is called a [fixed point](https://en.wikipedia.org/wiki/Fixed_point_(mathematics)).
This nomenclature, as well as the motivation behind this formulation, will become clear later on.

!!! info "Practice exercise"
	Convince yourself that any solution $$\bar{y}$$ to the equation
	$$0 = f(y)$$
	is also a solution to the equation
	$$y = g(y)$$,
	and vice versa.

The key idea is that when the fixed point $$\bar{y}$$ is __*input*__ into the function $$g(y)$$, the __*output*__ is the same quantity $$\bar{y}$$.
We then hope$$^\ddagger$$ that if we input a value __*close*__ to $$\bar{y}$$, the output will be a value that is __*even closer.*__
If this is true, then __*repeated substitution*__ will generate a sequence of values that approaches $$\bar{y}$$!
In mathematical terms, we wish to generate a sequence of approximate solutions
```math
\big\{ \, y_0, \, y_1, \, y_2, \, \ldots \, \big\}
```
where
```math
y_{n+1}
\, = \, g(y_n)
\qquad
\text{and}
\qquad
\lim_{n \rightarrow \infty} y_n
\, = \, \bar{y}
~.
```
This method is accordingly called [fixed point iteration](https://en.wikipedia.org/wiki/Fixed-point_iteration).

>  $$^\ddagger$$ rather, we will derive under what conditions repeated substitution will work

"""

# ╔═╡ 1b46f2a1-8ecb-47d4-9b94-d32dfdd359b0
md"""

### Attempt A (first attempt)

"""

# ╔═╡ 6caa8ec0-5d5d-4409-9db7-90f3aed8b4b6
gA(y) = y - f(y)

# ╔═╡ 547560c3-05dd-44ba-9374-e74d8d9e3d9b
yA0 = 1.5

# ╔═╡ 5f7548bf-dcf7-43a7-a35c-dae345751e2d
f(yA0) # want this value to be close to zero

# ╔═╡ aa38032b-f63b-4912-b6d6-1dae54487700
yA1 = gA(yA0)

# ╔═╡ 584a65ab-b575-4ac8-b40e-944909275b5b
f(yA1) # want this value to be CLOSER to zero

# ╔═╡ 58cc9fc2-f50a-4b04-98da-c7076d6a3164
yA2 = gA(yA1)

# ╔═╡ 51a3280b-9446-4e5f-b6b9-4ef51ad7f954
f(yA2) # want this value to be EVEN CLOSER to zero

# ╔═╡ 62211422-a046-4762-afce-855714e95ec8
md"""
!!! warning "Warning"
	The sequence $$\{ \, y^A_0, \, y^A_1, \, y^A_2 \, \}$$ is yielding __*progressively worse estimates*__ of $$\bar{y}$$!
	This means that __*the method is failing*__
"""

# ╔═╡ d4bef7c4-b8bc-4408-9dab-65370b7a8530
md"""

### Attempt B (second attempt)

The following trick is based on the error analysis calculation, which we will do next.
We note that if $$\bar{y}$$ satisfies
$$0 = 3 - y^2$$,
then it also satisfies
$$0 = (3 - y^2) / 2$$.
Let us then define a new function $$g^B(y)$$ as
```math
g^B(y)
\, = \, y
\, - \, \dfrac{1}{2} \, f(y)
\, = \, \dfrac{3}{2}
\, + \, y
\, - \, \dfrac{y^2}{2}
~.
```
Again,
$$\bar{y} = \sqrt{3}$$ is a solution to the equation
```math
y
\, = \, g^B(y)
~.
```
As before, we have reformulated the original problem (determine $$y$$ such that $$y^2 - 3 = 0$$) into a fixed-point problem.
Let us repeat our calculation.


"""

# ╔═╡ db5691bf-77da-4d75-9158-4a28f9a01bce
gB(y) = y - 0.5 * f(y)

# ╔═╡ 7a6ef620-7a10-454c-9b7c-15aa416b089f
yB0 = 1.5

# ╔═╡ 529db029-2b55-441e-990c-a04c2157c15e
f(yB0) # want this value to be close to zero

# ╔═╡ 07171d40-f399-4e2b-bec0-90761ee01ab7
yB1 = gB(yB0)

# ╔═╡ ede95cb9-67e3-4389-b732-819f1d09afb5
f(yB1) # want this value to be closer to zero

# ╔═╡ ebefd5e0-527c-4995-90d4-b8e148c7a7b1
yB2 = gB(yB1)

# ╔═╡ 43e14413-bb31-4f20-a22b-f8c2851715cc
f(yB2) # want this value to be even closer to zero!

# ╔═╡ 649991cf-b226-419a-9d58-54661214db0f
md"""

!!! question "Important Question"
	The only difference between the functions $$g(y)$$ and $$g^B(y)$$ is a factor of two.
	Why does one fixed point formulation work, while the other fails?
"""

# ╔═╡ f2e90575-23f8-425b-b2df-fe792a62ada5
md"""

### Error analysis

We seek to determine whether the sequence
```math
\big\{ \, y_0, \, y_1, \, y_2, \, \ldots \, \big\}
~,
\qquad
\text{where}
\qquad
y_{n+1}
\, = \, g(y_n)
~,
```
converges to the solution $$\bar{y}$$ of the equation
$$y = g(y)$$.
To this end, let us define $$r_n$$ as the residual for the approximate solution $$y_n$$:
```math
r_n
\, \equiv \, g(y_n)
\, - \, y_n
~.
```
In order for the series to converge, we require
```math
\lim_{n \rightarrow \infty} \big\lvert r_n \big\rvert
\, = \, 0
~.
```

> Our objective now is to express $$r_{n+1}$$ in terms of $$r_n$$.

As a starting point, we know
```math
r_n
\, = \, g(y_n)
\, - \, y_n
\qquad
\text{and}
\qquad
r_{n+1}
\, = \, g(y_{n+1})
\, - \, y_{n+1}
~.
```
Moreover, since
$$y_{n+1} = g(y_n)$$,
we can express $$r_n$$ as
```math
r_n
\, = \, y_{n+1}
\, - \, y_n
~.
```
__*Assuming $$r_n$$ is small,*__ i.e. our estimate $$y_n$$ is close to $$\bar{y}$$, we Taylor expand $$g(y_{n+1})$$ about $$y_n$$ and neglect higher-order terms as
```math
\begin{aligned}
g(y_{n+1})
\, &= \, g(y_n)
\, + \, (y_{n+1} - y_n) \, g'(y_n)
\, + \, \mathcal{O} (y_{n+1} - y_n)^2
\\[8pt]
\, &\approx \, g(y_n)
\, + \, (y_{n+1} - y_n) \, g'(y_n)
~.
\end{aligned}
```
Since $$g(y_n) = y_{n+1}$$, we find
```math
g(y_{n+1})
\, \approx \, y_{n+1}
\, + \, (y_{n+1} - y_n) \, g'(y_n)
~.
```
Substituting this result into the expression for $$r_{n+1}$$ yields
```math
r_{n+1}
\, \approx \, g'(y_n) \, r_n
~.
```
Accordingly, we expect our series to converge if
```math
\big\lvert g'(y_n) \big\rvert
\, < \, 1
\qquad
\text{for all}
\quad n
~,
```
as this guarantees
```math
\big\lvert r_{n+1} \big\rvert
\, < \, \big\lvert r_n \big\rvert
~.
```


"""

# ╔═╡ c7923532-c6d8-4666-878c-1bd2eb96e624
md"""

### Pros and cons

- simple implementation ✔
 

- convergence criterion more complicated when $$\boldsymbol{y}$$ and $$\boldsymbol{f}$$ are vectors ✖
- often slow to converge ✖

"""

# ╔═╡ 5a9839f4-c23c-497c-ae15-cbc5041f975f
md"""
## Newton–Raphson iteration

The final solution method we will consider is Newton–Raphson iteration, also sometimes called [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).
The method is founded on the Taylor series approximation __to linear order,__ and has a nice graphical representation in one dimension.
The basic idea will be demonstrated with our goal of calculatin $$\sqrt{3}$$ by finding the zero of the function
```math
f(y)
\, = \, y^2
\, - \, 3
~.
```
To use this method, we will need to calculate the derivative of the function under consideration:
```math
\dfrac{\mathrm{d} f}{\mathrm{d} y}
\, = \, 2 y
~.
```
Now, let us make an initial guess. We know that $$f(1) < 0$$ and $$f(2) > 0$$, so $$\sqrt{3}$$ lies somewhere between 1 and 2.
Let us choose as our initial guess
```math
y_0
\, = \, 1
~.
```

"""

# ╔═╡ 28a82140-a2a2-4c73-83cf-a7990ad25fde
df(y) = 2 * y

# ╔═╡ ac82c457-854e-43db-801a-128bddb85c82
y_vals = 0.0:0.1:3.0;

# ╔═╡ 8e7bc7f7-4e67-4aea-8023-9c9aac309c18
y0_NR = 1.0

# ╔═╡ b4d12ea4-b4a5-41f4-82f6-8a95b959bfc3
Plots.plot(
	[y_vals, [y0_NR]],
	[f.(y_vals), [f(y0_NR)]],
	label = ["f(y)" "f(y_0)"],
	xlabel = "y",
	xlims = (0.0, 3.0),
	ylims = (-3.0, 6.0),
	line = ([3 0], [:darkorange3 :deepskyblue3]),
	marker = ([:none :circle], [0 8], [:darkorange3 :deepskyblue3]),
	framestyle = :origin,
)

# ╔═╡ a34ddf42-a785-41c1-9ddc-cc466b80a84f
md"""
Clearly, this choice is not the correct answer.
According to the Newton–Raphson method, we calculate the next estimate in the sequence by taking a __linear approximation__ about the previous guess:
```math
f(y)
\, \approx \, \textcolor{Cyan}{\underbrace{f(y_0)}_{-2}}
\, + \, (y - y_0) \, \textcolor{Purple}{\underbrace{f'(y_0)}_{2}}
~.
```

"""

# ╔═╡ 2bdaed71-8697-456d-b5bc-0b04bbecdb64
f(y0_NR)

# ╔═╡ dfd90d98-aea1-4717-a958-164b4cd9478e
df(y0_NR)

# ╔═╡ fe180aa1-8de8-48c5-9c65-2eb6fd19702a
Plots.plot!(
	y_vals,
	(y_vals .- 1.0) * 2 .- 2,
	line = (2, :purple, :dash),
	label = "f0_approx",
)

# ╔═╡ b8df916b-7239-4143-8f9a-240ae5c29869
md"""
The linear approximation about $$y_0$$ intersects the horizontal axis at the point $$y = 2$$.
The Newton–Raphson method says to use this value as the next guess in the sequence:
```math
y_1
\, = \, 2
~.
```

!!! question "Important Question"
	How do we calculate where the linear approximation intersects the horizontal axis?

!!! info "Crucial Observation"
	The linear approximation is __zero__ at the intersection point!
	We thus find
	```math
	\textcolor{Red}{0}
	\, = \, f(y_0)
	\, + \, (\textcolor{Red}{y_1} - y_0) \, f'(y_0)
	~.
	```
With some rearrangement, we find
```math
y_1
\, = \, y_0
\, - \, \dfrac{f(y_0)}{f'(y_0)}
\, = \, 1
\, - \, \dfrac{-2}{2}
\, = \, 2
~.
```

"""

# ╔═╡ 942132be-9004-4a4e-be32-74f89b25ee9d
y1_NR = 2.0

# ╔═╡ ce377a5b-9b58-4e1c-91fe-caa86f24b12f
Plots.plot!(
	[y1_NR],
	[f(y1_NR)],
	label = "f(y_1)",
	line = (0, :red),
	marker = (:diamond, 8, :red),
)

# ╔═╡ 75779bdb-acd5-411b-a13b-7603cad61c70
Plots.plot!(
	y_vals,
	(y_vals .- 2.0) * 4 .+ 1,
	line = (2, :red, :dashdot),
	label = "f1_approx",
)

# ╔═╡ 21330e9e-b0c0-4e3c-b90e-501b3752aac4
md"""
> You'll notice that we are now much closer to the true solution, where the graph of $$f(y)$$ crosses the horizontal axis.
> We will continue iterating until we reach a user-specified threshold!
> Let's zoom in closer

"""

# ╔═╡ 2bce9785-8583-46ad-8f1a-ea8d3d04ebe3
Plots.plot(
	[y_vals, [y1_NR], y_vals],
	[f.(y_vals), [f(y1_NR)], (y_vals .- 2.0) * 4 .+ 1],
	label = ["f(y)" "f(y_1)" "f1_approx"],
	xlabel = "y",
	ylabel = "f(y)",
	xlims = (1.6, 2.1),
	ylims = (-0.7, 1.6),
	line = ([3 0 2], [:darkorange3 :red :red], [:solid :solid :dashdot]),
	marker = ([:none :diamond :none], [0 8 0], [:darkorange3 :red :red]),
	framestyle = :zerolines,
)

# ╔═╡ 48d9d911-bf39-4495-ad15-31455087eae8
md"""
We now see how after only two Newton–Raphson iterations, we are already close to the true solution!
The general procedure for the Newton–Raphson method is to generate a sequence
```math
\big\{
y_0, \, y_1, \, y_2, \, \ldots
\big\}
```
by choosing an initial value $$y_0$$ and then determining successive values according to
```math
y_{n+1}
\, = \, y_n
\, - \, \dfrac{f(y_n)}{f'(y_n)}
~.
```
We stop the iteration at a chosen tolerance, i.e. when
```math
f(y_n)
\, < \, \epsilon_{\text{tol}}
~.
```
"""

# ╔═╡ 554653bd-0b1a-4655-92a4-a9ccfee288c4
md"""

### Pros and cons

- extremely fast convergence ✔


- requires a good initial guess ✖
- requires calculation of the derivative ✖

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Plots = "~1.40.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "cd0152cf4ca63e41e1e8e7cfd4831fc9f5ed5873"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e16271d212accd09d52ee0ae98956b8a05c4b626"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "17.0.6+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═1f42dd9e-a82f-498d-a6b2-73d1a05603fd
# ╟─8268c060-58b4-11ef-0d04-dd88ba82f850
# ╠═a349ab6c-6d43-4850-9848-21d1d2158882
# ╠═1ffa3210-1442-499c-a066-1a9bb174f6a4
# ╠═e37983bb-04e5-4d01-bb15-9700774014c1
# ╟─aca8795a-a924-4718-b266-aec2a648a049
# ╠═9b937cc8-f70b-476f-b194-18bd5ca7fc91
# ╠═151957bb-cb42-43e8-9d9f-240b1776eaaf
# ╠═2df83b27-7dcf-4293-b390-c8434bb2a8f9
# ╟─ac0791ee-059b-4f03-b984-78bc90a2667a
# ╠═3028cddf-d42f-49af-9753-57f0dd1b55eb
# ╟─fe19ad8f-6275-4248-9e70-e3a8d6cf0bde
# ╠═9fb354dd-e0b3-4a49-b81d-af91e679c799
# ╟─de00f330-a035-4adb-bc88-3cd1c55edf75
# ╠═44c90555-ce2d-435b-93ad-ae58bee105c4
# ╟─72a4f34e-47dd-4be3-aaf7-ace5bd47f166
# ╠═24ef0a11-86cd-44d9-8695-4fa5480b9b0e
# ╠═5dc003d7-bea8-4744-9f0f-6a16151dcd9e
# ╠═f474de65-1416-4e7f-9855-5bc920e7314a
# ╟─f036f742-e5c7-45d2-91a6-3ccb0fe890be
# ╠═98540be7-1c08-4d39-9ecc-1d9bac36c913
# ╟─f7571334-7267-46e3-9418-1fb5aea32f54
# ╟─f2f1a24d-7e4c-4bda-ade3-101c455d2e51
# ╟─1b46f2a1-8ecb-47d4-9b94-d32dfdd359b0
# ╠═6caa8ec0-5d5d-4409-9db7-90f3aed8b4b6
# ╠═547560c3-05dd-44ba-9374-e74d8d9e3d9b
# ╠═5f7548bf-dcf7-43a7-a35c-dae345751e2d
# ╠═aa38032b-f63b-4912-b6d6-1dae54487700
# ╠═584a65ab-b575-4ac8-b40e-944909275b5b
# ╠═58cc9fc2-f50a-4b04-98da-c7076d6a3164
# ╠═51a3280b-9446-4e5f-b6b9-4ef51ad7f954
# ╟─62211422-a046-4762-afce-855714e95ec8
# ╟─d4bef7c4-b8bc-4408-9dab-65370b7a8530
# ╠═db5691bf-77da-4d75-9158-4a28f9a01bce
# ╠═7a6ef620-7a10-454c-9b7c-15aa416b089f
# ╠═529db029-2b55-441e-990c-a04c2157c15e
# ╠═07171d40-f399-4e2b-bec0-90761ee01ab7
# ╠═ede95cb9-67e3-4389-b732-819f1d09afb5
# ╠═ebefd5e0-527c-4995-90d4-b8e148c7a7b1
# ╠═43e14413-bb31-4f20-a22b-f8c2851715cc
# ╟─649991cf-b226-419a-9d58-54661214db0f
# ╟─f2e90575-23f8-425b-b2df-fe792a62ada5
# ╟─c7923532-c6d8-4666-878c-1bd2eb96e624
# ╟─5a9839f4-c23c-497c-ae15-cbc5041f975f
# ╠═28a82140-a2a2-4c73-83cf-a7990ad25fde
# ╠═ac82c457-854e-43db-801a-128bddb85c82
# ╠═8e7bc7f7-4e67-4aea-8023-9c9aac309c18
# ╠═b4d12ea4-b4a5-41f4-82f6-8a95b959bfc3
# ╟─a34ddf42-a785-41c1-9ddc-cc466b80a84f
# ╠═2bdaed71-8697-456d-b5bc-0b04bbecdb64
# ╠═dfd90d98-aea1-4717-a958-164b4cd9478e
# ╠═fe180aa1-8de8-48c5-9c65-2eb6fd19702a
# ╟─b8df916b-7239-4143-8f9a-240ae5c29869
# ╠═942132be-9004-4a4e-be32-74f89b25ee9d
# ╠═ce377a5b-9b58-4e1c-91fe-caa86f24b12f
# ╠═75779bdb-acd5-411b-a13b-7603cad61c70
# ╟─21330e9e-b0c0-4e3c-b90e-501b3752aac4
# ╠═2bce9785-8583-46ad-8f1a-ea8d3d04ebe3
# ╟─48d9d911-bf39-4495-ad15-31455087eae8
# ╟─554653bd-0b1a-4655-92a4-a9ccfee288c4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
