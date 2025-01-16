>[IMPORTANT]
>This code is only grammally check so it is not guaranteed physically correctly enough! Check it before you use it.

## Introduction
This is a code of a plasma simulation using the particle-in-cell (PIC) method. PIC is a method of simulating plasma with giant particles and discretizing fields into meshes[^1]. This code aims to simulate collisionless plasma with an ion-electron PIC code based on Matlab. The code is initialized as a hydrogen electron code, the user can change ion mass and charge by changing `mi` and `e_i`. This code allows injection of new particles into the chamber, user can define `particle_number_setted ` to limit maxima numbers of giant particles. Note that number of real particles contained in one giant particle is also determined by this number, set it carefully or it may influence result of the simulation.
[^1]: Dawson, J.M. (1983). "Particle simulation of plasmas". Reviews of Modern Physics. 55 (2): 403–447. Bibcode:1983RvMP...55..403D. doi:10.1103/RevModPhys.55.403.
