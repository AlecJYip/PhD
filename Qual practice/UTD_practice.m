clc;
clear;
close all;

f = 9e5;
c = 3e8;
lambda = c/f;
beta = 2*pi/lambda;

%beta = 1;

phi_p = 75*pi/180;
phi = 250*pi/180;
rho = 3*lambda;

% Derived constants
phi_p_m = phi + phi_p;

L = 3*lambda;

n = 2;

gamma_p = 90*pi/180;

[D_par, D_per] = UTD(L, phi, phi_p, n, beta,  gamma_p);