# about simulation code
## equation of motion
Equation of motion using in this simulation is as below.

![equation_of_motion](../images/equation_of_motion.png)

## parameters
Parameters using in this simulation is as below

|parameter|value|
|:--:|:--:|
|initial altitude [km]|300|
|magnetic latitude [°]|66|
|electric field|0|
|magnetic field|detail is as below|
|energy of particle(electron) [keV]|1|
|mirror force|with／ignore|
|initial pitch angle [°]|0-90|
|width of time [s]|3.07 ×  10^(-8)|
|time step [times]|2.0 ×  10^7|
|collision probability|detail is as below|

### about magnetic field
Magnetic field is approximated with dipole field.

![magnetic_field](../images/magnetic_field.png)

![r_0](../images/surface_of_earth.png): Earth radius

![B_0](../images/magnetic_field_of_r0.png): magnitude of magnetic field at Earth's surface



##### [reference] _BasicSpacePlasmaPhysics, p32, eq3.1_
equation of dipole field is as below.

![equation_of_dipole_field](../images/dipole_magnetic_field.png)


##### [reference] _Katoh and Omura, JGR 2006_
In this simulation, divergence of magnetic field is set to 0

![divergenceB](../images/divergence_of_magnetic_field.png)


### about collision probability
Collicion probability is refered to graph of collision frequency in _BasicSpacePlasmaPhysics,p.66_

mean of number of collision during cyclotron period, λ is led as below approximated with the graph.
(75km is standard)

![mean_of_collision](../images/mean_of_collision.png)


When collision occurs according to exponential distribution,
collision probability, P per 1 step in this simulation is lead as below.

![collision_probability](../images/collision_probability.png)



### about normalization
In this simulation, parameters are normalized.

Unit of parameters is as below

| parameter | unit |
|:----------:|:---:|
|velocity|c^(-1)|
|time|Ω^(-1)|
|distance|cΩ^(-1)|
|mass|m\_e^(-1)|
|charge|q\_e^(-1)|

c: speed of light,
Ω: cyclotron period,
m\_e: mass of electron,
q\_e: charge of electron

In this simulation, Ω is value of it at Earth's surface.
