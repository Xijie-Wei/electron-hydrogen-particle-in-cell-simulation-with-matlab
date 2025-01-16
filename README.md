>[!IMPORTANT]
>This code is only grammally checked so it is not guaranteed physically correctly enough! Check it before you use it.

## Introduction
This is a code of a plasma simulation using the particle-in-cell (PIC) method. PIC is a method of simulating plasma with giant particles and discretizing fields into meshes[^1]. This code aims to simulate collisionless plasma with an ion-electron PIC code based on Matlab. The code is initialized as a hydrogen electron code, the user can change ion mass and charge by changing `mi` and `e_i`. This code allows injection of new particles into the chamber, user can define `particle_number_setted ` to limit maxima numbers of giant particles. Note that number of real particles contained in one giant particle is also determined by this number, set it carefully or it may influence the result of the simulation.

Here is the code for the initialization of the particles. `particle_number_e` and `particle_number_i` are used to set the initial number of electrons and ions, respectively. `particle_e` and `particle_i` are used to set the initial position of electrons and ions. `v_drift` is used to set an initial drift velocity of the particles, can be used for simulate plasma source. 
```
%Intialze the particle
random_number = rand(particle_number_e,3);
random_number_v = rand(particle_number_e,3)-0.5;
particle_e = [(0*random_number(:,1))*x_total,(0.3*random_number(:,2)+0.35)*y_total,(0.1*random_number(:,3)+0.45)*z_total];
particle_v_e = zeros(particle_number_e,3);
for i = 1:particle_number_e
    particle_v_e(i,1) = v_drift+v0_e*random_number_v(i,1)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
    particle_v_e(i,2) = v0_e*random_number_v(i,2)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
    particle_v_e(i,3) = v0_e*random_number_v(i,3)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
end

random_number = rand(particle_number_i,3);
random_number_v = rand(particle_number_i,3)-0.5;
particle_i = [(0*random_number(:,1))*x_total,(0.3*random_number(:,2)+0.35)*y_total,(0.1*random_number(:,3)+0.45)*z_total];
particle_v_i = zeros(particle_number_i,3);
for i = 1:particle_number_i
    particle_v_i(i,1) = v_drift+v0_i*random_number_v(i,1)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
    particle_v_i(i,2) = v0_i*random_number_v(i,2)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
    particle_v_i(i,3) = v0_i*random_number_v(i,3)/sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
end
```
[^1]: Dawson, J.M. (1983). "Particle simulation of plasmas". Reviews of Modern Physics. 55 (2): 403–447. Bibcode:1983RvMP...55..403D. doi:10.1103/RevModPhys.55.403.
