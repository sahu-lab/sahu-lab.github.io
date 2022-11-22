---
layout: post
title:  "Nosé–Hoover thermostat: Theory"
author: "Amaresh Sahu"
date:   2022-11-22 11:00:00 -0000
---

The Nosé--Hoover thermostat refers to a class of deterministic algorithms for [molecular dynamics][Molecular_dynamics]{:target="_blank"} (MD) simulations at constant temperature.
In such simulations, the system of interest interacts with a [heat bath][Thermal_reservoir], and the dynamics of both the system and the reservoir are calculated over time.
The thermostat is named after [Shūichi Nosé][Shuichi_Nose]{:target="_blank"}, who in 1984 proposed the first method for simulations of the coupled system and heat bath,[^ref_nose_mp]<sup>,</sup>[^ref_nose_jcp] and [William G. Hoover][William_G_Hoover]{:target="_blank"}---who in the following year improved Nosé's method and presented the coupled dynamics in the form commonly used today.[^ref_hoover_pra]
A remarkable feature of the Nosé--Hoover equations is that only a single fictitious particle is used to describe the heat bath.
Consequently, the Nosé--Hoover thermostat is computationally inexpensive, and is widely used in constant temperature molecular dynamics simulations.


# Why do we need a thermostat?

A "traditional" molecular dynamics simulation will sample the [microcanonical ensemble][Microcanonical_ensemble]{:target="_blank"}, in which the total energy is fixed, rather than the [canonical ensemble][Canonical_ensemble]{:target="_blank"}, in which the temperature is fixed.
To see why, we consider a simulation involving $$ N $$ particles occupying a volume $$ V $$, where the position and momentum of the $$ \alpha^{\textrm{th}} $$ particle are denoted $$ \boldsymbol{q}^\alpha $$ and $$ \boldsymbol{p}^\alpha $$, respectively.
The [(micro)state][Microstate]{:target="_blank"} of a classical system is captured by the set of all particle positions and momenta, for which we introduce the shorthand

$$
	\boldsymbol{\Gamma}
	\, := \, \big\{\mkern1mu
		\boldsymbol{q}^1, \,
		\boldsymbol{q}^2, \,
		\ldots \,
		\boldsymbol{q}^N; \
		\boldsymbol{p}^1, \,
		\boldsymbol{p}^2, \,
		\ldots \,
		\boldsymbol{p}^N
	\mkern1mu\big\}
	~.
$$

If the system is [isolated][Isolated_system]{:target="_blank"} and governed by a Hamiltonian $$ \mathcal{H} = \mathcal{H} (\boldsymbol{\Gamma}) $$, then the positions and momenta evolve in time according to [Hamilton's equations][Hamiltonian_mechanics]{:target="_blank"}---written as

$$
	\boldsymbol{\dot{q}}^\alpha
	\, = \, \dfrac{\partial \mathcal{H}}{\partial \boldsymbol{p}^\alpha}
	\qquad\quad
	\text{and}
	\qquad\quad
	\boldsymbol{\dot{p}}^\alpha
	\, = \, - \, \dfrac{\partial \mathcal{H}}{\partial \boldsymbol{q}^\alpha}
	~,
$$

where a 'dot' accent denotes a time derivative.
In this way, the system traces out a trajectory $$ \boldsymbol{\Gamma} (t) $$ in [phase space][Phase_space].
By taking the total time derivative of the Hamiltonian and then substituting the equations of motion, we find

$$
\dot{\mathcal{H}}
	\ = \, \sum_{\alpha = 1}^N \bigg(
		\boldsymbol{\dot{p}}^\alpha \boldsymbol{\cdot} \mkern1mu \dfrac{\partial \mathcal{H}}{\partial \boldsymbol{p}^\alpha}
		\, + \, \boldsymbol{\dot{q}}^\alpha \boldsymbol{\cdot} \mkern1mu \dfrac{\partial \mathcal{H}}{\partial \boldsymbol{q}^\alpha}
	\bigg)
	\, = \, \sum_{\alpha = 1}^N \bigg(
		\boldsymbol{\dot{p}}^\alpha \boldsymbol{\cdot} \boldsymbol{\dot{q}}^\alpha
		\, - \, \boldsymbol{\dot{q}}^\alpha \boldsymbol{\cdot} \boldsymbol{\dot{p}}^\alpha
	\bigg)
	\, = \ 0
	~.
$$

Accordingly, the value of the Hamiltonian---namely, the [total (internal) energy][Internal_energy]{:target="_blank"} $$ E $$---is a constant of motion.
We have thus shown that **an MD simulation will sample the microcanonical ensemble, in which $$ N $$, $$ V $$, and $$ E $$ are constant, if particles evolve in time according to Hamilton's equations.**

The development of microcanonical MD simulations marks a significant achievement in statistical physics.
However, the results of such simulations often cannot directly be compared to experimental measurements---in which the temperature $$ T $$, rather than the internal energy $$ E $$, is the relevant control parameter.
There were consequently many efforts to formulate algorithms for constant temperature MD simulations, one of which is the Nosé--Hoover thermostat.


# The development by Nosé

For an MD simulation to sample the canonical ensemble, the system of interest must be perturbed such that its total energy does not remain constant. Nosé's great insight was to construct an extended, fictitious system involving one additional degree of freedom that represented a heat bath.
**If the fictitious system is constructed appropriately, then microcanonical simulations of the fictitious system lead to the real system sampling the canonical ensemble.**
In this fashion, one can simulate Hamilton's equations in the *extended* system and guarantee the *real* system is at the desired temperature.
In what follows, we present the main ideas of Nosé's developments, albeit in a slightly different form from the original formulation.[^ref_nose_mp]<sup>,</sup>[^ref_nose_jcp]


We begin by constructing a fictitious system, where all associated quantities are denoted with a 'hat' accent.
At present, we do not know what these quantities represent, and we will determine their physical meaning subsequently.
It is most natural to begin with the fictitious [Lagrangian][Lagrangian_mechanics]{:target="_blank"} $$ \hat{\mathcal{L}} $$, which is a function of
$$ \{ \boldsymbol{\hat{q}}^\alpha \} $$,
$$ \{ \mathrm{d} \boldsymbol{\hat{q}}^\alpha / \mathrm{d} \hat{t} \} $$,
$$ \hat{s} $$,
and
$$ \mathrm{d} \hat{s} / \mathrm{d} \hat{t} $$,
and is given by

$$
	\hat{\mathcal{L}}
	\, = \, \sum_{\alpha = 1}^N \dfrac{m \hat{s}^2}{2} \, \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}} \boldsymbol{\cdot} \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\, - \, U \big(
		\{\boldsymbol{\hat{q}}^\alpha\}
	\big)
	\, + \, \dfrac{Q}{2} \bigg(
		\dfrac{\mathrm{d} \hat{s}}{\mathrm{d} \hat{t}}
	\bigg)^2
	- \, \big(
		3 N + 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln \hat{s}
	~.
$$

Here $$ k_{\mathrm{B}} $$ is [Boltzmann's constant][Boltzmann_constant], $$ T_{\mathrm{ext}} $$ is the temperature we seek to impose externally, and $$ U $$ is the [potential energy][Potential_energy]{:target="_blank"} of the system.
The quantity $$ \hat{s} $$ is dimensionless, implying the corresponding inertial term $$ Q $$ has dimensions of energy $$ \cdot $$ time<sup>2</sup>.
However, since $$ \hat{s} $$ is thought of as relating to a heat bath, $$ Q $$ is often called the *thermal mass* or *thermal inertia*.
With the Lagrangian, the [conjugate momenta][Generalized_momentum]{:target="_blank"} are determined according to the standard procedure:

$$
	\boldsymbol{\hat{p}}^\alpha
	\, = \, \dfrac{\partial \hat{\mathcal{L}}}{\partial ( \mathrm{d} \boldsymbol{\hat{q}}^\alpha / \mathrm{d} \hat{t} \mkern1mu) }
	\, = \, m \hat{s}^2 \, \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\qquad\quad
	\text{and}
	\qquad\quad
	\hat{p}^{}_{\! s}
	\, = \, \dfrac{\partial \hat{\mathcal{L}}}{\partial ( \mathrm{d} \hat{s} / \mathrm{d} \hat{t} \mkern1mu) }
	\, = \, Q \, \dfrac{\mathrm{d} \hat{s}}{\mathrm{d} t}
	~.
$$

We can then construct the corresponding fictitious Hamiltonian as the Legendre transform of the Lagrangian---expressed as

$$
	\hat{\mathcal{H}} \big(
		\{ \boldsymbol{\hat{q}}^\alpha \}, \,
		\{ \boldsymbol{\hat{p}}^\alpha \}, \,
		\hat{s}, \,
		\hat{p}^{}_{\! s}
	\big)
	\, = \, \sum_{\alpha = 1}^N \, \boldsymbol{\hat{p}}^\alpha \boldsymbol{\cdot} \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\, + \, \hat{p}^{}_{\! s} \, \dfrac{\mathrm{d} \hat{s}}{\mathrm{d} \hat{t}}
	\, - \, \hat{\mathcal{L}}
	~,
$$

from which we find

$$
	\hat{\mathcal{H}}
	\, = \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{\hat{p}}^\alpha \boldsymbol{\cdot} \boldsymbol{\hat{p}}^\alpha}{2 m \hat{s}^2}
	\, + \, U \big(
		\{\boldsymbol{\hat{q}}^\alpha\}
	\big)
	\, + \, \dfrac{p^2_{\! s}}{2 Q}
	\, + \, \big(
		3 N + 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln \hat{s}
	~.
$$

At this point, two tasks remain.
We seek to interpret all of the fictitious quantities and relate them to our physical system, and then show that microcanonical simulations of the extended system (as governed by $$ \hat{\mathcal{H}} $$) cause the physical system to sample the canonical distribution.


## The physical interpretation of fictitious quantities

It is immediately clear from the fictitious Lagrangian that we need to interpret $$ \hat{s} $$, as if $$ \hat{s} = 1 $$ then $$ \hat{\mathcal{L}} $$ is the usual physical Lagrangian governing an $$ N $$-particle system at constant energy.
To this end, we examine how the fictitious Lagrangian is modified under the choice of scaling

$$
	\hat{s}
	\, \rightarrow \, b \mkern1mu \hat{s}
	~,
	\qquad
	\hat{t}
	\, \rightarrow \, b^m \mkern1mu \hat{t}
	~,
	\qquad
	\text{and}
	\qquad
	\hat{\boldsymbol{q}}^\alpha
	\, \rightarrow \, b^n \mkern1mu \hat{\boldsymbol{q}}^\alpha
	~,
$$

where the exponents $$ m $$ and $$ n $$ will be chosen such that the Lagrangian is invariant to the value of the dimensionless scale factor $$ b $$.
Under such a change of variables, the fictitious Lagrangian transforms to

$$
	\hat{\mathcal{L}}
	\, \rightarrow \, b^{2 - 2m + 2n} \mkern1mu \sum_{\alpha = 1}^N \dfrac{m \hat{s}^2}{2} \, \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}} \boldsymbol{\cdot} \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\, - \, U \big(
		\{b^n \mkern1mu \boldsymbol{\hat{q}}^\alpha\}
	\big)
$$

$$
	\hspace{20pt} + \, b^{2 - 2m} \mkern1mu \dfrac{Q}{2} \bigg(
		\dfrac{\mathrm{d} \hat{s}}{\mathrm{d} \hat{t}}
	\bigg)^2
	- \, \big(
		3 N + 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln (b \mkern1mu \hat{s})
	~.
	\hspace{-20pt}
$$

If we choose $$ m = 1 $$ and $$ n = 0 $$, then the Lagrangian is only shifted by the constant factor $$ - (3N + 1) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln b $$---which does not affect the dynamics.
Since $$ \boldsymbol{\hat{q}}^\alpha $$ is unaffected by the scaling, we interpret it as the real position:
$$ \boldsymbol{\hat{q}}^\alpha = \boldsymbol{q}^\alpha $$.
We also observe that
$$ \hat{s} \rightarrow b \mkern1mu \hat{s} $$
and
$$ \hat{t} \rightarrow b \mkern1mu \hat{t} $$,
for which the ratio $$ \hat{t} / \hat{s} $$ is invariant to the scaling.
Since $$ \hat{t} $$ corresponds to the real time when $$ \hat{s} = 1 $$, we interpret $$ \hat{t} $$ as a *scaled time*, with $$ \hat{s} $$ the (continuously changing) *value* of the scaling.
Since MD simulations of the fictitious system will generate $$ \hat{s} $$ as a function of $$ \hat{t} $$, we can convert between a fictitious time $$ \hat{\tau} $$ and real time $$ \tau $$ according to

$$
	\mathrm{d} t
	\, = \, \dfrac{\mathrm{d} \hat{t}}{\hat{s}}
	\qquad
	\text{and}
	\qquad
	\tau
	\, = \, \int_0^{\hat{\tau}} \dfrac{\mathrm{d} \hat{t}}{\hat{s} (\hat{t})}
	~.
$$



## The sampling of the canonical ensemble

Our final task is to show how solving Hamilton's equations in the extended, fictitious system sets the temperature to a constant value in the real system.
We do so by examining the [partition function][Partition_function]{:target="_blank"}, which describes how often a system can be found in a particular microstate.

In the microcanonical ensemble, the partition function $$ W (E) $$ is simply the number of microstates in which the system has total energy $$ E $$---with each state being equally probable.
For the fictitious system, we have

$$
	\hat{W} (\hat{E})
	\, = \, \dfrac{E_0}{N! \, h^{3N+1}} \,
	\int \mathrm{d} \hat{p}^{}_{\! s}
	\int \mathrm{d} \hat{s} \,
	\prod_{\alpha = 1}^N \bigg(
		\int \mathrm{d} \boldsymbol{\hat{p}}^\alpha
		\int \mathrm{d} \boldsymbol{\hat{q}}^\alpha
	\bigg) \,
	\bigg\{
		\delta \Big(
			\hat{\mathcal{H}} \big(
				\{ \boldsymbol{\hat{q}}^\alpha \}, \,
				\{ \boldsymbol{\hat{p}}^\alpha \}, \,
				\hat{s}, \,
				\hat{p}^{}_{\! s}
			\big)
			\, - \, \hat{E}
		\mkern1mu \Big)
	\bigg\}
	~,
$$

where $$ \delta (\ldots) $$ is the [Dirac $$ \delta $$-function][Dirac_delta_function]{:target="_blank"}, $$ h $$ is [Planck's constant][Planck_constant]{:target="_blank"}, and $$ E_0 $$ is an arbitrary constant with units of energy.
Note the factor of $$ E_0 / h^{3N + 1} $$ ensures $$ \hat{W} $$ is dimensionless, and the factor of $$ N! $$ accounts for the $$ N $$ particles being indistinguishable.
At this point, Nosé recognized it was useful to *recast* the partition function in terms of the real positions and momenta of the particles, rather than their fictitious counterparts.
We already noted
$$ \boldsymbol{\hat{q}}^\alpha = \boldsymbol{q}^\alpha $$,
and now recognize

$$
	\boldsymbol{\hat{p}}^\alpha
	\, = \, m \hat{s}^2 \, \dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\, = \, m \hat{s} \, \dfrac{\mathrm{d} \boldsymbol{q}^\alpha}{\mathrm{d} t}
	~,
	\qquad
	\text{for which}
	\qquad
	\boldsymbol{p}^\alpha
	\, = \, \dfrac{\boldsymbol{\hat{p}}^\alpha}{\hat{s} \,}
	\, = \, m \mkern1mu \boldsymbol{\dot{q}}^\alpha
$$

is the physical momentum.
The fictitious Hamiltonian is then expressed in terms of the physical microstate $$ \boldsymbol{\Gamma} $$ as

$$
	\hat{\mathcal{H}} \big(
		\boldsymbol{\Gamma}, \,
		\hat{s}, \,
		\hat{p}^{}_{\! s}
	\big)
	\, = \, \mathcal{H} ( \boldsymbol{\Gamma} )
	\, + \, \dfrac{p^2_{\! s}}{2 Q}
	\, + \, \big(
		3 N + 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln \hat{s}
	~,
$$

where

$$
	\mathcal{H} ( \boldsymbol{\Gamma} )
	\, = \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{p}^\alpha \boldsymbol{\cdot} \boldsymbol{p}^\alpha}{2 m}
	\, + \, U \big(
		\{\boldsymbol{q}^\alpha\}
	\big)
$$

is the physical Hamiltonian.
Moreover, we recognize
$$ \mathrm{d} \boldsymbol{\hat{p}}^\alpha = \hat{s}^3 \mathrm{d} \boldsymbol{p}^\alpha $$
in a three-dimensional system, and express the partition function as

$$
	\hat{W} ( \hat{E} )
	\, = \, \dfrac{E_0}{N! \, h^{3N+1}}
	\int \! \mathrm{d} \hat{p}^{}_{\! s}
	\int \! \mathrm{d} \hat{s}
	\int \! \mathrm{d} \boldsymbol{\Gamma} \,
		\hat{s}^{3 N} \, \delta \bigg(
			\mathcal{H} ( \boldsymbol{\Gamma} )
			\, + \, \dfrac{p^2_{\! s}}{2 Q}
			\, + \, \big(
				3 N + 1
			\big) \mkern1mu k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \mkern1mu \ln \hat{s}
			\, - \, \hat{E}
		\mkern1mu \bigg)
	~.
$$

Note that here we introduced the shorthand
$$
	\mathrm{d} \boldsymbol{\Gamma} := \prod_{\alpha = 1}^N (
		\int \mathrm{d} \boldsymbol{q}^\alpha
		\int \mathrm{d} \boldsymbol{p}^\alpha
	)
$$.
Importantly, the argument of the $$ \delta $$-function now contains only a single term involving $$ \hat{s} $$.
By applying the $$ \delta $$-function identity

$$
	\delta \big( \hat{g} (\hat{s}) \big)
	\, = \, \dfrac{\delta ( \hat{s} - \hat{s}_0 )}{\hat{g}{\mkern1mu}'(\hat{s}_0)}
$$

for any function $$ \hat{g} (\hat{s}) $$ with
$$
	g (\hat{s}_0) = 0
$$,
and introducing
$$
	\beta_{\mathrm{ext}}
	:= ( k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} )^{-1}
$$
as the [inverse temperature][Thermodynamic_beta]{:target="_blank"}, we find

$$
	\delta \bigg(
		\mathcal{H} ( \boldsymbol{\Gamma} )
		\, + \, \dfrac{p^2_{\! s}}{2 Q}
		\, + \, \big(
			3 N + 1
		\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \ln \hat{s}
		\, - \, \hat{E}
	\mkern1mu \bigg)
	\, = \, \dfrac{\hat{s}_0 \mkern1mu \beta_{\mathrm{ext}}}{3N + 1}
	\, \delta \big(
		\hat{s}
		\, - \, \hat{s}_0
	\big)
	~,
$$

where

$$
	\hat{s}_0
	\, = \, \exp \bigg\{
		\dfrac{\beta_{\mathrm{ext}}}{3N + 1} \, \bigg(
			\hat{E}
			\, - \, \hat{\mathcal{H}} ( \boldsymbol{\Gamma} )
			\, - \, \dfrac{p^2_{\! s}}{2 Q}
		\bigg)
	\bigg\}
	~.
$$

The integration over $$ \hat{s} $$ is now easily carried out.
After some rearrangement, we obtain

$$
	\hat{W} (\hat{E})
	\, = \, \dfrac{
		\beta^{}_{\mathrm{ext}} \mkern1mu E_0 \,
		\exp ( \beta^{}_{\mathrm{ext}} \mkern1mu \hat{E} )
	} {(3N + 1) \, N! \, h^{3N+1}} \,
	\int \mathrm{d} \hat{p}^{}_{\! s} \exp \bigg\{
		- \dfrac{\beta^{}_{\mathrm{ext}} \, p_{\! s}^2}{(2 Q)}
	\bigg\}
	\int \mathrm{d} \boldsymbol{\Gamma} \,
	\exp \Big\{
		- \beta^{}_{\mathrm{ext}} \mkern1mu \mathcal{H} ( \boldsymbol{\Gamma} )
	\Big\}
	~.
$$

The [Gaussian integral][Gaussian_integral]{:target="_blank"} over $$ \hat{p}^{}_{\! s} $$ can be calculated exactly, with which the fictitious microcanonical partition function is expressed as

$$
	\hat{W} (\hat{E})
	\, = \, \dfrac{
		E_0 \, \sqrt{ 2 \pi Q \beta^{}_{\mathrm{ext}} \,} \,
		\exp (
			\beta^{}_{\mathrm{ext}} \mkern1mu \hat{E}
		)
	}{(3N + 1) \, h} \,
	Q^{}_{ N V T_{\mathrm{ext}}}
	~.
$$

Here,

$$
	Q^{}_{ N V T_{\mathrm{ext}}}
	\, := \, \dfrac{1}{N! \, h^{3N}} \,
	\int \mathrm{d} \boldsymbol{\Gamma} \,
	\exp \Big\{
		- \beta^{}_{\mathrm{ext}} \mkern1mu \mathcal{H} ( \boldsymbol{\Gamma} )
	\Big\}
$$

is the [*canonical partition function*][Canonical_partition_function]{:target="_blank"} of a system of $$ N $$ particles in a volume $$ V $$, maintained at a constant temperature $$ T_{\mathrm{ext}} $$.


## Significance

The two previous equations contain a crucial result from the analysis of Nosé, and it is useful to further highlight their importance.
Let us imagine simulating the fictitious system at a constant energy $$ \hat{E} $$: each microstate with
$$
	\hat{\mathcal{H}} (\boldsymbol{\hat{q}}^\alpha, \boldsymbol{\hat{p}}^\alpha, \hat{s}, \hat{p}^{}_{\! s})
	\, = \, \hat{E}
$$
is equally likely.
However, if we then *only examine* the physical degrees of freedom, we find the probability $$ \mathcal{P} $$ of observing a physical microstate $$ \boldsymbol{\Gamma} $$ satisfies

$$
	\mathcal{P} (\boldsymbol{\Gamma})
	\, \sim \, \exp \Big\{
		- \beta_{\mathrm{ext}} \mkern1mu \mathcal{H} (\boldsymbol{\Gamma})
	\Big\}
	~,
$$

which is *exactly* the probability distribution of a system at the desired temperature $$ T_{\mathrm{ext}} $$.
**We have thus shown that when the fictitious system is in the microcanonical ensemble, the real system samples the canonical ensemble at the desired temperature $$ T_{\mathrm{ext}} $$.**


## Practical considerations

The above result concludes our presentation of Nose's theoretical developments.
In practice, standard MD methods to simulate the fictitious system in the microcanonical ensemble would generate

$$
	\boldsymbol{\hat{q}}^\alpha (\hat{t}), \,
	\qquad
	\boldsymbol{\hat{p}}^\alpha (\hat{t}), \,
	\qquad
	\hat{s} (\hat{t}), \,
	\qquad
	\text{and}
	\qquad\quad
	\hat{p}^{}_{\! s} (\hat{t})
$$

at a set of discrete fictitious times
$$
	\hat{t}
	\in \{
		\hat{t}_0, \,
		\hat{t}_1, \,
		\hat{t}_2, \,
		\ldots, \,
		\hat{t}_T
	\}
$$.
Generally, the time interval $$ \Delta \hat{t} $$ in simulations is constant, for which
$$
	\hat{t}_j
	= j \Delta \hat{t}
$$.
From the scale factor $$ \hat{s} $$ at discrete $$ \hat{t} $$, we approximate the real time according to

$$
	t_i
	\, = \, \int_0^{t_i} \mathrm{d} t
	\, = \, \int_0^{\hat{t}_i} \dfrac{\mathrm{d} \hat{t}}{\hat{s} (\hat{t})}
	\, \approx \, \sum_{j = 0}^i \dfrac{w_j \, \Delta \hat{t}}{\hat{s} (\hat{t}_j)}
	~,
$$

where the $$ \{ w_j \} $$ are the weights of some [quadrature formula][Gaussian_quadrature]{:target="_blank"}---for example, the [trapezoidal rule][Trapezoidal_rule]{:target="_blank"}.
The real positions and momenta are then given by

$$
	\boldsymbol{q}^\alpha (t_i)
	\, = \, \boldsymbol{\hat{q}}^\alpha (\hat{t}_i)
	\qquad
	\text{and}
	\qquad
	\boldsymbol{p}^\alpha (t_i)
	\, = \, \dfrac{\boldsymbol{\hat{p}}^\alpha (\hat{t}_i)}{ \hat{s} (\hat{t}_i) }
	~,
$$

with which all other observables can be calculated.
Importantly, the calculated $$ \boldsymbol{\Gamma} (t_i) $$ appropriately sample the canonical ensemble at the specified temperature $$ T_{\mathrm{ext}} $$.


# The extension by Hoover

Nosé's equations mark a significant advancement in constant temperature MD simulations.
However, the equations are slightly inconvenient to use in practice because of the need to calculate $$ t $$ from $$ \hat{t} $$, and the resultant non-uniform (real) time intervals $$ \Delta t $$.
In the year following Nosé's seminal work, Hoover recognized it was possible to simulate Nosé's equations in terms of the *real* (rather than fictitious) time---and do away with any time scaling altogether.[^ref_hoover_pra]
In what follows, we present the main results of Hoover's theoretical developments.

We begin with the equations of motion for the extended, fictional system.
By substituting the fictional Hamiltonian into Hamilton's equations of motion, we find

$$
	\dfrac{\mathrm{d} \boldsymbol{\hat{q}}^\alpha}{\mathrm{d} \hat{t}}
	\, = \, \dfrac{\boldsymbol{\hat{p}}^\alpha}{m \hat{s}^2}
	~,
	\qquad\qquad
	\dfrac{\mathrm{d} \boldsymbol{\hat{p}}^\alpha}{\mathrm{d} \hat{t}}
	\, = \, - \, \dfrac{\partial U}{\partial \boldsymbol{\hat{q}}^\alpha}
	~,
	\qquad\qquad
	\dfrac{\mathrm{d} \hat{s}}{\mathrm{d} \hat{t}}
	\, = \, \dfrac{\hat{p}^{}_{\! s}}{Q}
	~,
$$

and

$$
	\dfrac{\mathrm{d} \hat{p}^{}_{\! s}}{\mathrm{d} \hat{t}}
	\, = \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{\hat{p}}^\alpha \boldsymbol{\cdot} \boldsymbol{\hat{p}}^\alpha}{m \hat{s}^3}
	\, - \, \dfrac{(3N+1) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}}}{\hat{s}}
	~.
$$

At this point, Hoover recognized the utility in expressing these equations entirely in terms of real quantities.
To this end, we introduce the following fictitious quantities that depend on the *real* time (note that we drop the 'hat' accent):

$$
	\eta
	\, := \, \ln \hat{s}
	\qquad\quad
	\text{and}
	\qquad\quad
	p^{}_{\! \eta}
	\, := \, \hat{p}^{}_{\! s}
	~.
$$

By recognizing

$$
	\dfrac{\mathrm{d} \hat{s}}{\mathrm{d} \hat{t}}
	\, = \, \dfrac{1}{\hat{s}} \, \dfrac{\mathrm{d} \hat{s}}{\mathrm{d} t}
	\, = \, \dfrac{\mathrm{d} \ln \hat{s}}{\mathrm{d} t}
	\, = \, \dfrac{\mathrm{d} \eta}{\mathrm{d} t}
	\, = \, \dot{\eta}
$$

and defining

$$
	\boldsymbol{f}^\alpha
	\, := \, - \, \dfrac{\partial U}{\partial \boldsymbol{q}^\alpha}
$$

as the force on the $$ \alpha^{\mathrm{th}} $$ particle, we express the governing equations in terms of the real time $$ t $$ as

$$
	\boldsymbol{\dot{q}}^\alpha
	\, = \, \dfrac{\boldsymbol{p}^\alpha}{m}
	~,
	\qquad\quad
	\boldsymbol{\dot{p}}^\alpha
	\, = \, \boldsymbol{f}^\alpha
	\, - \, \bigg(
		\dfrac{p^{}_{\! \eta}}{Q}
	\bigg) \, \boldsymbol{p}^\alpha
	~,
	\qquad\quad
	\dot{\eta}
	\, = \, \dfrac{p^{}_{\! \eta}}{Q}
	~,
$$

and

$$
	\dot{p}^{}_{\! \eta}
	\, = \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{p}^\alpha \boldsymbol{\cdot} \boldsymbol{p}^\alpha}{m}
	\, - \, (3N+1) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}}
	~.
$$

In this manner, the dynamics of the fictitious system are simulated as a function of the real time $$ t $$, rather than the fictitious time $$ \hat{t} $$.
Consequently, the conversion between $$ \hat{t} $$ and $$ t $$ is no longer required, and an MD simulation generates

$$
	\boldsymbol{q}^\alpha (t)
	~,
	\qquad
	\boldsymbol{q}^\alpha (t)
	~,
	\qquad
	\eta (t)
	~,
	\qquad
	\text{and}
	\qquad
	p^{}_{\! \eta} (t)
$$

at a set of discrete and uniformly spaced *real* times
$$
	t_j
	= j \mkern1mu \Delta t
$$.
The practical consequences of this formulation are significant, and the above equations---referred to as the Nosé--Hoover equations---are generally used in numerical implementations.


## Additional insights

With the Nosé--Hoover equations expressed in terms of the real time $$ t $$, we can draw additional insights about the dynamics of the extended system.
First, the quantity $$ (p^{}_{\! \eta} / Q) $$ acts as a drag coefficient: each physical particle feels an additional force $$ - (p^{}_{\! \eta} / Q) \boldsymbol{p}^\alpha $$ due to the heat bath.
In this case, however, the drag coefficient dynamically evolves in time and can be either positive or negative.
Moreover, recall from the [equipartition thoerem][Equipartition_theorem]{:target="_blank"} that

$$
	k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}}
	\, = \, \dfrac{1}{N_{\mathrm{DOFs}}} \, \bigg\langle
		\sum_{\alpha = 1}^N \dfrac{\boldsymbol{p}^\alpha \boldsymbol{\cdot} \boldsymbol{p}^\alpha}{m}
	\bigg\rangle
	~,
$$

where $$ N_{\mathrm{DOFs}} $$ is the number of physical [degrees of freedom][Degrees_of_freedom]{:target="_blank"} and $$ \langle \ldots \rangle $$ denotes an [ensemble average][Ensemble_average]{:target="_blank"}.
Consequently, we define an "*instantaneous temperature*" $$ T(t) $$ according to

$$
	k_{\mathrm{B}} \mkern1mu T(t)
	\, := \, \dfrac{1}{N_{\mathrm{DOFs}}} \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{p}^\alpha \boldsymbol{\cdot} \boldsymbol{p}^\alpha}{m}
	~,
$$

for which the dynamical equation governing $$ p^{}_{\! \eta} $$ can be written as

$$
	\dot{p}^{}_{\! \eta}
	\, = \, N_{\mathrm{DOFs}} \, k_{\mathrm{B}} \mkern1mu T(t)
	\, - \, \big(
		3 N
		+ 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}}
	~.
$$

If $$ N_{\mathrm{DOFs}} $$ were to equal $$ 3N + 1 $$, then there would be a simple physical interpretation for the drag force, which affects all particle momenta:

- it "cools down" a system that is too "hot," for which $$ T(t) > T_{\mathrm{ext}} $$
- it "heats up" a system that is too "cold," for which $$ T(t) < T_{\mathrm{ext}} $$

As it turns out,
$$
	N_{\mathrm{DOFs}}
	\ne 3 N + 1
$$.
However, many numerical implementations of the Nosé--Hoover thermostat in fact *replace* the factor of $$ (3N + 1) $$ with $$ N_{\mathrm{DOFs}} $$ due to the physical motivation described above.[^ref_allen_tildesley]
This subtlety will be discussed in more detail in a future post.

We close by noting that the Nosé equations can be derived from the fictitious Hamiltonian $$ \mathcal{H} (\boldsymbol{\hat{q}}^\alpha, \boldsymbol{\hat{p}}^\alpha, \hat{s}, \hat{p}^s) $$.
However, the same is **not** true for the Nosé--Hoover equations, in which the fundamental variables are the real (rather than fictitious) particle positions and momenta.
Nonetheless, the *value* of the fictitious Hamiltonian remains a constant of motion in the extended system, and can be expressed in terms of real quantities.
We thus define the quantity

$$
	E_{\mathrm{NH}}
	\, := \, \sum_{\alpha = 1}^N \dfrac{\boldsymbol{p}^\alpha \boldsymbol{\cdot} \boldsymbol{p}^\alpha}{2 m}
	\, + \, U \big(
		\{\boldsymbol{q}^\alpha\}
	\big)
	\, + \, \dfrac{p^2_{\! \eta}}{2 Q}
	\, + \, \big(
		3 N + 1
	\big) \, k_{\mathrm{B}} \mkern1mu T_{\mathrm{ext}} \, \eta
	~,
$$

which has the dimensions of an energy and is conserved by the Nosé--Hoover equations.


# Numerical details

Since the Nosé--Hoover equations are non-Hamiltonian, they cannot be solved numerically using the usual integration methods from constant energy MD simulations.
The numerical simulation of the Nosé--Hoover equations will be discussed in a future post.
In the meantime, we recommend the excellent textbook by Allen and Tildesley.[^ref_allen_tildesley]


# References

[^ref_nose_mp]: S. Nosé. <a href="https://doi.org/10.1080/00268978400101201" target="_blank">A molecular dynamics method for simulations in the canonical ensemble</a>. *Mol. Phys.* **52**, 255--268 (1984)

[^ref_nose_jcp]: S. Nosé. <a href="https://doi.org/10.1063/1.447334" target="_blank">A unified formulation of the constant temperature molecular dynamics methods</a>. *J. Chem. Phys.* **81**, 511--519 (1984)

[^ref_hoover_pra]: W.G. Hoover. <a href="https://doi.org/10.1103/PhysRevA.31.1695" target="_blank">Canonical Dynamics: Equilibrium phase-space distributions</a>. *Phys. Rev. A* **31**, 1695--1697 (1985)

[^ref_allen_tildesley]: M.P. Allen and D.J. Tildesley. <a href="https://doi.org/10.1093/oso/9780198803195.001.0001" target="_blank">Computer Simulation of Liquids</a>, Chapter 3.8.2, 2<sup>nd</sup> Edition. New York: Oxford University Press, 2017. See also the <a href="https://global.oup.com/booksites/content/9780198803195/" target="_blank">Companion Website</a> and <a href="https://github.com/Allen-Tildesley/examples" target="_blank">Github Repository</a>.


<!-- Wikipedia links -->

[Boltzmann_constant]: https://en.wikipedia.org/wiki/Boltzmann_constant
[Canonical_ensemble]: https://en.wikipedia.org/wiki/Canonical_ensemble
[Canonical_partition_function]: https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)#Canonical_partition_function
[Degrees_of_freedom]: https://en.wikipedia.org/wiki/Degrees_of_freedom_(physics_and_chemistry)
[Dirac_delta_function]: https://en.wikipedia.org/wiki/Dirac_delta_function
[Equipartition_theorem]: https://en.wikipedia.org/wiki/Equipartition_theorem
[Ensemble_average]: https://en.wikipedia.org/wiki/Ensemble_(mathematical_physics)#Ensemble_average
[Gaussian_integral]: https://en.wikipedia.org/wiki/Gaussian_integral
[Gaussian_quadrature]: https://en.wikipedia.org/wiki/Gaussian_quadrature
[Generalized_momentum]: https://en.wikipedia.org/wiki/Momentum#Generalized
[Internal_energy]: https://en.wikipedia.org/wiki/Internal_energy
[Isolated_system]: https://en.wikipedia.org/wiki/Isolated_system
[Hamiltonian_mechanics]: https://en.wikipedia.org/wiki/Hamiltonian_mechanics
[Lagrangian_mechanics]: https://en.wikipedia.org/wiki/Lagrangian_mechanics
[Microcanonical_ensemble]: https://en.wikipedia.org/wiki/Microcanonical_ensemble
[Microstate]: https://en.wikipedia.org/wiki/Microstate_(statistical_mechanics)
[Molecular_dynamics]: https://en.wikipedia.org/wiki/Molecular_dynamics
[Nose_Hoover]: https://en.wikipedia.org/wiki/Nos%C3%A9%E2%80%93Hoover_thermostat
[Partition_function]: https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)
[Phase_space]: https://en.wikipedia.org/wiki/Phase_space
[Planck_constant]: https://en.wikipedia.org/wiki/Planck_constant
[Potential_energy]: https://en.wikipedia.org/wiki/Potential_energy
[Shuichi_Nose]: https://en.wikipedia.org/wiki/Shuichi_Nos%C3%A9
[Thermal_reservoir]: https://en.wikipedia.org/wiki/Thermal_reservoir
[Thermodynamic_beta]: https://en.wikipedia.org/wiki/Thermodynamic_beta
[Trapezoidal_rule]: https://en.wikipedia.org/wiki/Trapezoidal_rule
[William_G_Hoover]: https://en.wikipedia.org/wiki/William_G._Hoover



