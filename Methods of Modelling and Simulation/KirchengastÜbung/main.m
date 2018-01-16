% Übung Kirchengast
% Felix Sonnleitner, 01430166
clc;
clear all;
close all;

Re = 6371;  %km
z = 650;    %km
theta0 = 0;
phi0 = 0;   %30, 60, 90, 89,9
N0 = 315;
Hn = 7;     %km
r0 = Re;

span = Re:0.1:Re+z