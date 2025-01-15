%%% This code is used for 3-D particle in cell plasma simulation
clear;
close all;
%%% Constants will be defined as below
global miu0 ep0
miu0 = 1.25663706127e-6 % N⋅A−2
ep0 = 8.8541878188e-12 % F⋅m−1
me = 9.1093837e-31 % kg
mi = 1.67e-27 % kg
e0 = 1.602176634e-19 % C
kb = 1.380649e-23 %J⋅K−1

%User set constant
T0 = 1000 %K
ne = 1e14 %m-3
global grid_size
grid_size = [40,40,2] %number of grid in each direction
total_step = 500 %total step of the simulation
global particle_number_i particle_number_e
particle_number_setted =1e5% Setted number of particle used in this simulation
particle_number_i = 200 %Number of ion used in this simulation
particle_number_e = 200 %Number of electron used in this simulation
v_drift = 1e5;
max_injection = 200 %Maxima number of particle injected in every step

%Calculated constant
global debye_length
debye_length = sqrt(ep0*kb*T0/(ne*e0^2)) %m
plasma_frequency = sqrt(ne*e0^2/(me*ep0)) %rad/s
global t_unit
t_unit = 0.02*2*pi / plasma_frequency %s time unit of the simutaion
x_total = grid_size(1) * debye_length %Physical x length simulated
y_total = grid_size(2) * debye_length %Physical x length simulated
z_total = grid_size(3) * debye_length %Physical x length simulated
total_number = x_total*y_total*z_total*ne %Total electron number
global number_represent
number_represent = total_number / particle_number_setted %Number of electrons a particle presents
e_p = -e0 * number_represent %Charge of a electron
m_p_e = me * number_represent %Mass of a electron
e_i = e0 * number_represent %Charge of a ion
m_p_i = mi * number_represent %Mass of a ion
v0_e = sqrt(kb*T0*3/(me*number_represent))
v0_i = sqrt(kb*T0*3/(mi*number_represent))

function B = maxwell_get_B(B0,E)
%%% This function is used to calculate magnatic field
%%% The function returns a 3-D array contains information of B
%%% The equation used in this function is dB/dt = - curl(E)
global t_unit
[curlx,curly,curlz,cav] = curl(E.x,E.y,E.z);
B.x = B0.x - curlx * t_unit;
B.y = B0.y - curly * t_unit;
B.z = B0.z - curlz * t_unit;
end

function E = maxwell_get_E(E0,B,J)
%%% This function is used to calculate electical field
%%% The function returns a 3-D array contains information of E
%%% The equation used in this function is dE/dt = (curl(B)/miu0 - J)/ep0
global miu0 ep0 t_unit
[curlx,curly,curlz,cav] = curl(B.x,B.y,B.z);
E.x = E0.x + t_unit * (curlx/miu0 - J.x) / ep0;
E.y = E0.y + t_unit * (curly/miu0 - J.y) / ep0;
E.z = E0.z + t_unit * (curlz/miu0 - J.z) / ep0;
end

function n = number_density(number,position)
%%% This function is used to calculate number density in the grid of
%%% electon
%%% The function returns a 3-D array contains information of number
%%% density in grid
global  grid_size number_represent debye_length
n = zeros(grid_size);
for i = 1 : number
    [posi,weight] = weighting(position(i,:));
    for j = 1:length(weight)
        n(posi(j,1),posi(j,2),posi(j,3)) = n(posi(j,1),posi(j,2),posi(j,3)) + weight(j)*number_represent/debye_length^3;
    end
end
end

function [position,weight] = weighting(r)
%%% This function is used to calculated contribution of a particle to the
%%% number density, field and J
%%% This function returns a array with position vector and weighting
global debye_length grid_size

posi(:) = (floor(r(:)/debye_length));
d_p(:) =   r(:)/debye_length - posi(:);

position(1,:) = [posi(1),posi(2),posi(3)+1];
weight(1) = (1-d_p(1))*(1-d_p(2))*d_p(3);
position(2,:) = [posi(1),posi(2)+1,posi(3)+1];
weight(2) = (1-d_p(1))*d_p(2)*d_p(3);
position(3,:) = [posi(1),posi(2)+1,posi(3)];
weight(3) = (1-d_p(1))*d_p(2)*(1-d_p(3));
position(4,:) = [posi(1),posi(2),posi(3)];
weight(4) = (1-d_p(1))*(1-d_p(2))*(1-d_p(3));

position(5,:) = [posi(1),posi(2),posi(3)+1];
weight(5) = d_p(1)*(1-d_p(2))*d_p(3);
position(6,:) = [posi(1)+1,posi(2)+1,posi(3)+1];
weight(6) = d_p(1)*d_p(2)*d_p(3);
position(7,:) = [posi(1)+1,posi(2)+1,posi(3)];
weight(7) = d_p(1)*d_p(2)*(1-d_p(3));
position(8,:) = [posi(1)+1,posi(2),posi(3)];
weight(8) = d_p(1)*(1-d_p(2))*(1-d_p(3));

for i = 1 : 8
    for j = 1:3
        if position(i,j) < 1
            position(i,j) = 1;
        end
        if position(i,j) > grid_size(j) - 1
            position(i,j) = grid_size(j) - 1;

        end
    end
end


% position(:) = (round(r(:)/debye_length));
% weight = 1;
% for i = 1:3
%     if position(i) < 1
%         position(i) = 1;
%     end
%     if position(i) > grid_size(i) - 1
%         position(i) = grid_size(i) - 1;
%
%     end
% end
end

%%% This is the part of this code

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

%Initalze the field
B.x = zeros(grid_size);
B.y = zeros(grid_size);
B.z = zeros(grid_size);
E.x = zeros(grid_size);
E.y = zeros(grid_size);
E.z = zeros(grid_size);
J.x = zeros(grid_size);
J.y = zeros(grid_size);
J.z = zeros(grid_size);

%User defined extiernal field
E_ext.x = zeros(grid_size);
E_ext.y = zeros(grid_size);
E_ext.z = zeros(grid_size);

B_ext.x = zeros(grid_size);
B_ext.y = 1e-3*ones(grid_size);
B_ext.z = zeros(grid_size);

f0 = 2450e6;
c = 3e9;

[x_mesh,y_mesh] = meshgrid(0:debye_length*1e3:(y_total-debye_length)*1e3,0:debye_length*1e3:(x_total-debye_length)*1e3);

%Simulation code
for t = 1 : total_step
    clc
    t
    t*t_unit
    y_total
    E = maxwell_get_E(E,B,J);
    % A laser bin injected
    for i = round(20:grid_size(1))
        for j = round(1:grid_size(3))
            E_ext.y(i,:,j) = (120)*sin((1:grid_size(2))/(c/f0)+ 2*pi*t*f0);
        end
    end

    J.x = zeros(grid_size);
    J.y = zeros(grid_size);
    J.z = zeros(grid_size);

    for j = 1 : particle_number_e
        [posi,weight]=weighting(particle_e(j,:));
        B_local = [0,0,0];
        E_local = [0,0,0];
        for k = 1 : length(weight)
            B_local(1) = B_local(1) + B.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(1) = E_local(1) + E.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            B_local(2) = B_local(2) + B.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(2) = E_local(2) + E.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            B_local(3) = B_local(3) + B.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(3) = E_local(3) + E.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
        end
        F(:) = e_p * (E_local(1,:)+cross(particle_v_e(j,:),B_local(1,:)));
        a(:) = F(:) / m_p_e;
        particle_v_e(j,:) = particle_v_e(j,:) + a(:).' * t_unit;
        particle_e(j,:) = particle_e(j,:) + particle_v_e(j,:) * t_unit;
    end

    for j = 1 : particle_number_i
        [posi,weight]=weighting(particle_i(j,:));
        B_local = [0,0,0];
        E_local = [0,0,0];
        for k = 1 : length(weight)
            B_local(1) = B_local(1) + B.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(1) = E_local(1) + E.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.x(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            B_local(2) = B_local(2) + B.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(2) = E_local(2) + E.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.y(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            B_local(3) = B_local(3) + B.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + B_ext.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
            E_local(3) = E_local(3) + E.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k) + E_ext.z(posi(k,1),posi(k,2),posi(k,3)) * weight(k);
        end
        F(:) = e_i * (E_local(1,:)+cross(particle_v_i(j,:),B_local(1,:)));
        a(:) = F(:) / m_p_i;
        particle_v_i(j,:) = particle_v_i(j,:) + a(:).' * t_unit;
        particle_i(j,:) = particle_i(j,:) + particle_v_i(j,:) * t_unit;
    end

    %kill particle out off boundary electron
    count = 0;
    particle_temp = zeros(particle_number_e,3);
    particle_v_temp = zeros(particle_number_e,3);
    for j = 1 : particle_number_e
        % z _axis reflected
        if particle_e(j,3) > z_total
            particle_e(j,3) = z_total - (particle_e(j,3) - z_total);
            particle_v_e(j,3) = - particle_v_e(j,3);
        end
        if particle_e(j,3) < 0
            particle_e(j,3) = - particle_e(j,3);
            particle_v_e(j,3) = - particle_v_e(j,3);
        end

        if not(particle_e(j,1) < 0 || particle_e(j,1) > x_total ||particle_e(j,2) < 0 || particle_e(j,2) > y_total ||particle_e(j,3) < 0 || particle_e(j,3) > z_total)
            particle_temp(count+1,:) = particle_e(j,:);
            particle_v_temp(count+1,:) = particle_v_e(j,:);
            count = count + 1;
        end
    end
    particle_e = particle_temp;
    particle_v_e = particle_v_temp;
    particle_number_e = count

     %kill particle out off boundary ion
    count = 0;
    particle_temp = zeros(particle_number_i,3);
    particle_v_temp = zeros(particle_number_i,3);
    for j = 1 : particle_number_i
        % z _axis reflected
        if particle_i(j,3) > z_total
            particle_i(j,3) = z_total - (particle_i(j,3) - z_total);
            particle_v_i(j,3) = - particle_v_i(j,3);
        end
        if particle_i(j,3) < 0
            particle_i(j,3) = - particle_i(j,3);
            particle_v_i(j,3) = - particle_v_i(j,3);
        end

        if not(particle_i(j,1) < 0 || particle_i(j,1) > x_total ||particle_i(j,2) < 0 || particle_i(j,2) > y_total ||particle_i(j,3) < 0 || particle_i(j,3) > z_total)
            particle_temp(count+1,:) = particle_i(j,:);
            particle_v_temp(count+1,:) = particle_v_i(j,:);
            count = count + 1;
        end
    end
    particle_i = particle_temp;
    particle_v_i = particle_v_temp;
    particle_number_i = count
    
    f=figure(1);
    f.Position=[100,100,800,300];
    subplot(1,2,1);
    d_n = number_density(particle_number_e,particle_e);
    den = sum(d_n,3);
    hold on
    pcolor(x_mesh,y_mesh,den,edgecolor='none');
    contour(x_mesh,y_mesh,den,5,'w');
    hold off
    xlabel("x / mm")
    ylabel("y / mm")
    title("Electron density")
    colorbar()
    subplot(1,2,2);
    d_n = number_density(particle_number_i,particle_i);
    den = sum(d_n,3);
    hold on
    pcolor(x_mesh,y_mesh,den,edgecolor='none');
    contour(x_mesh,y_mesh,den,5,'w');
    hold off
    xlabel("x / mm")
    ylabel("y / mm")
    title("Ion density")
    colorbar()

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

    for j = 1 : particle_number_i
        [posi,weight] = weighting(particle_i(j,:));
        for k = 1 : length(weight)
            J.x(posi(k,1),posi(k,2),posi(k,3)) = J.x(posi(k,1),posi(k,2),posi(k,3)) + e_i * particle_v_i(j,1) * weight(k);
            J.y(posi(k,1),posi(k,2),posi(k,3)) = J.y(posi(k,1),posi(k,2),posi(k,3)) + e_i * particle_v_i(j,2) * weight(k);
            J.z(posi(k,1),posi(k,2),posi(k,3)) = J.z(posi(k,1),posi(k,2),posi(k,3)) + e_i * particle_v_i(j,3) * weight(k);
        end
    end

    B = maxwell_get_B(B,E);
end
