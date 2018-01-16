% Übung Kirchengast
% Felix Sonnleitner, 01430166
clc;
clear all;
close all;

Re     = 6371;           %km
z      = 650;            %km
theta0 = 0;              %Grad
phi0   = [0 30 60 89.9]; %Grad
N0     = 315;            %N-units
Hn     = 7;              %km
r0     = Re;             %km

begin = r0;
stop  = Re + z;
step  = 0.1;

[theta_euler, r_euler] = f_euler(phi0, Hn, N0, theta0, begin, step, stop);
[theta_rk, r_rk]       = f_rungeKutta(phi0, Hn, N0, theta0, begin, step, stop);

hold on 
plot(r_euler,theta_euler)
plot(r_rk,theta_rk)
hold off