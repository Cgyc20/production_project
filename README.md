## Hybrid modelling the reaction diffusion system

This code will model the following system of production-degradation in one dimensional space:

$$
\emptyset \xrightarrow{k_1} A, \quad A \xrightarrow{k_2} \emptyset.
$$

Here we have rate of production equal to $k_1$ and rate of degradation equal to $k_2$. We implement this over both time and space over the given domain:

$$
\Omega = [0,L], \quad  t \in [0,T].
$$

We will be using the reaction diffusion master equation to model the system using the Gillespie stochastic simulation algorithm, this involves dividing the domain into $M$ total compartments, given by $C_1, \cdots, C_M$. 

Diffusion is modelled by considering a jump-rate, $d$, between neighbouring compartments. Species $A$ has a local copy within each corresponding compartment:

$$
A_1 \xrightleftharpoons[\ d\ ]{d} A_2 \xrightleftharpoons[\ d\ ]{d} \cdots \xrightleftharpoons[\ d\ ]{d} A_M.
$$

This system equates to the following partial differential equation, in which $D = dh^2$, $h$ is the compartment width. We let $u(x,t)$ be a continuous function denoting the concentration of species $A$ at point $x$ in space, and $t$ in time.

$$
\frac{\partial c}{\partial t} = D\frac{\partial^2 c}{\partial x^2} + k_1 - k_2c.
$$

## Method of lines

We will use finite difference methods to approximate the numerical solution. We will use the method of lines to create a mesh in two-dimensional space in the $(x,t)$ axes. We define the following:

$$
x_j = j\Delta x, \ t_n=n\Delta t \quad j=0,1,\cdots ,J \quad i = 0,1,\cdots, N
$$
$$
\Delta x = \frac{L}{J}, \quad \Delta t = \frac{T}{N}
$$

In this case we will let the numerical solution $u^i_j \approx u(x_j,t_n)$. How we solve this will be discussed in the following section, this is because we solve a slightly different PDE.
## The hybrid method

The hybrid system will model the reacting species using both the stochastic simulation algorithm and the PDE. The PDE will be used in regions of high particle count, while the stochastic algorithm will be applied to regions with low particle count. The hybrid method will convert ,ass between either particle representation. This adaptively occurs across each individual compartment. We let there be $M$ stochastic compartments and $J = KM+1$ PDE points (defined over a finer domain). 

The important feature of the hybrid method is that we let there be two types of mass within the simulation. The **Discrete mass** be given by the following: $D(C_j,t) \ \forall j \in [1,M]$. The discrete mass has a local copy for each stochastic compartment, and handled by the SSA.

We will define continuous mass, $u(x,t)$ to denote the concentration of the species $D$ in continuous time and space. This particular PDE is given by the following:

$$
\frac{\partial c}{\partial t} = D\frac{\partial^2 c}{\partial x^2} + \cancel{k_1} - k_2c.
$$

Here we remove the production term, as in the hybrid model - we want all production being governed by the stochastic simulation algorithm. We then let the conversion reactions occur between the Discrete and continuous mass depending on the total combined particle number within that region. The problem here is that the discrete mass $D(C_j,t)$ is an integer value, while $u




$$
\int_{\mathcal{C}_i} c(x,t) \, dx = M_i
$$

We can approximate $M_i$ to be the following:

$$
\int_{\mathcal{C}_i} c(x,t) \, dx = \sum_{j=ik}^{j=(i+1)k-1} c(x_j,t) \Delta x
$$


We can then approximate the total mass to be $A_i = D_i+ M_i$ at each compartment. 

Note we will exclude the production term in the PDE, only the stochastic regime will produce mass.

