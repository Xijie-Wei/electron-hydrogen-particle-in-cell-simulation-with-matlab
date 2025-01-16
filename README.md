[>Important]: This code is only grammally check so it is not guaranteed physically correctly enough! Check it before you use it.

## Introduction
This is a code of a plasma simulation using the particle-in-cell (PIC) method. PIC is a method of simulating plasma with giant particles and discretizing fields into meshes[^1]. This code aims to simulate collisionless plasma with an ion-electron PIC code based on Matlab. The code is initialized as a hydrogen electron code, the user can change ion mass and charge by changing mi and e_i

`
mi = 1.67e-27 % kg
`
for ion mass and
`e_i = e0 * number_represent %Charge of a ion
`
for ion charge

[^1]: Dawson, J.M. (1983). "Particle simulation of plasmas". Reviews of Modern Physics. 55 (2): 403–447. Bibcode:1983RvMP...55..403D. doi:10.1103/RevModPhys.55.403.
