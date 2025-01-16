>[!IMPORTANT]
>This code is only grammally checked so it is not guaranteed physically correctly enough! Check it before you use it.

## Introduction
This is a code of a plasma simulation using the particle-in-cell (PIC) method. PIC is a method of simulating plasma with giant particles and discretizing fields into meshes[^1]. This code aims to simulate collisionless plasma with an ion-electron PIC code based on Matlab. 
## Code details for change set-ups of the simulation
Here are some important variables determine the set-up of the simulation
1. The code is initialized as a hydrogen electron code, the user can change ion mass and charge by changing `mi` and `e_i`.
2. This code allows injection of new particles into the chamber, user can set `particle_number_setted ` to limit the maxima numbers of giant particles and `max_injection` for the maxima number of giant particles injected each time step. Note that number of real particles contained in one giant particle is also determined by this number, set it carefully or it may influence the result of the simulation.
3. `T0` and `ne` are used to set the initial temperature and number density of the particles, in K and $m^{-3}$ respectively.
4. `grid_size` is a 3-D array containing number of grips in x,y,z axis, respectively.
5. `total_step` is the number of steps simulated

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
Here is the code for the injection of the new particles, following the similar procedure of the initialization of the particles
```
    %inject new particle electron
    inject_n = min(max_injection,particle_number_setted-particle_number_e)
    random_number = rand(inject_n,3);
    particle_e(particle_number_e+1:particle_number_e+inject_n,:) =...
        [(0*random_number(:,1))*x_total,(0.3*random_number(:,2)+0.35)*y_total,(0.1*random_number(:,3)+0.45)*z_total];
    random_number_v = rand(inject_n,3)-0.5;
    for i = 1 : inject_n
        leng = sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
        particle_v_e(i+particle_number_e,1) = v_drift+v0_e*random_number_v(i,1)/leng;
        particle_v_e(i+particle_number_e,2) = v0_e*random_number_v(i,2)/leng;
        particle_v_e(i+particle_number_e,3) = v0_e*random_number_v(i,3)/leng;
    end
    particle_number_e = particle_number_e + inject_n;

    %inject new particle ion
    inject_n = min(max_injection,particle_number_setted-particle_number_i)
    random_number = rand(inject_n,3);
    particle_i(particle_number_i+1:particle_number_i+inject_n,:) =...
        [(0*random_number(:,1))*x_total,(0.3*random_number(:,2)+0.35)*y_total,(0.1*random_number(:,3)+0.45)*z_total];
    random_number_v = rand(inject_n,3)-0.5;
    for i = 1 : inject_n
        leng = sqrt(random_number_v(i,1)^2+random_number_v(i,2)^2+random_number_v(i,3)^2);
        particle_v_i(i+particle_number_i,1) = v_drift+v0_i*random_number_v(i,1)/leng;
        particle_v_i(i+particle_number_i,2) = v0_i*random_number_v(i,2)/leng;
        particle_v_i(i+particle_number_i,3) = v0_i*random_number_v(i,3)/leng;
    end
    particle_number_i = particle_number_i + inject_n;

    %calculate J
    for j = 1 : particle_number_e
        [posi,weight] = weighting(particle_e(j,:));
        for k = 1 : length(weight)
            J.x(posi(k,1),posi(k,2),posi(k,3)) = J.x(posi(k,1),posi(k,2),posi(k,3)) + e_p * particle_v_e(j,1) * weight(k);
            J.y(posi(k,1),posi(k,2),posi(k,3)) = J.y(posi(k,1),posi(k,2),posi(k,3)) + e_p * particle_v_e(j,2) * weight(k);
            J.z(posi(k,1),posi(k,2),posi(k,3)) = J.z(posi(k,1),posi(k,2),posi(k,3)) + e_p * particle_v_e(j,3) * weight(k);
        end
    end
```
In the part of killing particles out of the boundary, user can apply codes to reflect particles or produce a continued plasma tube, here use the electron particle in z axis as an example. For reflecting particles hit the boundary.
```
% z _axis reflected
        if particle_e(j,3) > z_total
            particle_e(j,3) = z_total - (particle_e(j,3) - z_total);
            particle_v_e(j,3) = - particle_v_e(j,3);
        end
        if particle_e(j,3) < 0
            particle_e(j,3) = - particle_e(j,3);
            particle_v_e(j,3) = - particle_v_e(j,3);
        end
```
For continued plasma tube.
```
        if particle_e(j,3) > z_total
            particle_e(j,3) = particle_e(j,3) - z_total;
        end
        if particle_e(j,3) < 0
            particle_e(j,3) = particle_e(j,3) + z_total;
        end

```
For x axis and y axis, change `particle_e(j,3)` to `particle_e(j,1)` or `particle_e(j,2)`, `z_total` to `x_total` or `y_total`, respectively. Similar procedure for `particle_i()` to change boundary conditions for ion.
One important thing determines the simulation setup is external field, in this code, it is set a microwave is injected into the plasma.
```
f0 = 2450e6;
c = 3e9;
```
```
    for i = round(20:grid_size(1))
        for j = round(1:grid_size(3))
            E_ext.y(i,:,j) = (120)*sin((1:grid_size(2))/(c/f0)+ 2*pi*t*f0);
        end
    end
```
For other types of external fields, change `E_ext` and `B_ext`. e.g. `E_ext.x` determines x weight of external electrical field.
[^1]: Dawson, J.M. (1983). "Particle simulation of plasmas". Reviews of Modern Physics. 55 (2): 403–447. Bibcode:1983RvMP...55..403D. doi:10.1103/RevModPhys.55.403.
