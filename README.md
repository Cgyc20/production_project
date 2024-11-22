# Hybrid modelling the reaction diffusion system


Welcome to this Python project! This models a reaction diffusion system in one dimensional space. This project uses [Poetry](https://python-poetry.org/) as a dependency manager


## Prerequisites

Before running the project, ensure you have the following installed:

### 1. Python 3.12
Download and install Python 3.12 from the official [Python website](https://www.python.org/downloads/).

### 2. Poetry
Install Poetry, the dependency manager used in this project:


### 3. Clone the repository

```bash 
git clone https://github.com/Cgyc20/production_project.git
cd production_project
``` 

### 4. Install package dependencies
Run the following command to install the required dependencies in a virtual environment:

```bash 
poetry install
``` 

## Running the project

### Input parameters 



The project requires a file named `parameters_input.dat` to provide input configurations.


### Running scripts

To execute the main project scripts, use the following commands:


#### To run the main script:
```bash 
poetry run main
``` 

#### To animate the saved data (only after running the main script at least once):
```bash 
poetry run animate
``` 




# The Theory
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

Here we remove the production term, as in the hybrid model - we want all production being governed by the stochastic simulation algorithm. We then let the conversion reactions occur between the Discrete and continuous mass depending on the total combined particle number within that region. The problem here is that the discrete mass $D(C_j,t)$ is an integer value, while $u(x_i,t_j)$ is defined over a finer domain, and hence in order to convert mass we must workout the total mass of the continuous function over that particular compartment. For this we just do the following:

$$
\int_{C_j} u(x,t) \, dx = U_j
$$

We integrate the PDE over the spatial compartment to generate a total mass value of $U_j$. In this code we use the left hand rule to integrate across the PDE. Given we calculate this value, then we can calculate the total mass within the compartment as being the summation of either methods:

$$
A_j = D_j + U_j
$$

### Mass conversion

We control the conversion of mass between either regime using the following two reactions:

$$
D \xrightarrow{\gamma} U, \ \text{Conversion from Discrete to continuous},
$$
$$
U \xrightarrow{\gamma} D \ \text{Conversion from continuous to discrete}.
$$

Note that conversion to continuous only occurs when $A_i \geq \text{threshold}$ within that compartment, and the opposite for conversion to discrete. 

In order to conserve mass then we must convert this mass correctly. 

Given we have $D_j \xrightarrow{\gamma} U_j$ in compartment $j$, then we must have that a particles worth of mass is correctly added to the continuous regime. This involves adding $1/h$ particles worth of mass for each point in the corresponding domain of the PDE solution. 

