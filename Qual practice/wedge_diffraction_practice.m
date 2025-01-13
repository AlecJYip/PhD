clc;
clear;
close all;

% Constants
f = 10e9;
c = 3e8;
lambda = c/f;
beta = 2*pi/lambda;
phi_p = 75*pi/180;
phi = 250*pi/180;
rho = 3*lambda;

% Derived constants
phi_p_m = phi + phi_p;
alpha = 1 + cos(phi_p_m);
const_1 = -exp(1i*pi/4)*sqrt(2/(pi*alpha));
const_2 = exp(1i*beta*rho*cos(phi_p_m))*cos(phi_p_m/2);

% Integration limits and step size
upper_t = 1000;
lower_t = sqrt(alpha * beta * rho);
d_t = 0.0001;
% Vectorized calculation

%t = lower_t:d_t:upper_t; % Generate the vector of t values
t = lower_t:d_t:upper_t;
sume=sum(exp(-1i * t.^2)*d_t); % Vectorized summation


% Result
v_b = const_1 * const_2 * sume;


%% other

n = 2;

term_1 = exp(-1i*(beta*rho+pi/4))/sqrt(2*pi*beta*rho);

term_2 = (1/n)*sin(pi/n)/(cos(pi/n)-cos(phi_p_m/n));

v_b_wedge = term_1*term_2;


%% poli

term_1 = 2*exp(1i*pi/4)/(n*sqrt(pi));

term_2 = sin(pi/n)/(cos(pi/n)-cos(phi_p_m/n));

term_3 = abs(cos(phi_p_m/2))*exp(1i*beta*rho*cos(phi_p_m));

t = lower_t:d_t:upper_t;

v_b_poli=term_1*term_2*term_3*sum(exp(-1i * t.^2)*d_t); % Vectorized summation


