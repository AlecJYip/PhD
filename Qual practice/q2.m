clc;
clear;
close all;

c = 3e8;

f = 1e9;

lambda = c/f;

R = 10*lambda;

d_theta = pi/18000;

F = 0;

for theta = d_theta:d_theta:pi

  fun = func(theta);

  diff_theta = sin(theta)*d_theta;

  F = F + fun*diff_theta;

end

F = F*2*pi;

theta_end = pi;

theta_start = 0;


G_t = zeros(1,length(theta_start:d_theta:theta_end));

counter = 1;

for theta = theta_start:d_theta:theta_end

G_t(counter) = func(theta);

counter = counter + 1;

end


G_t2 = zeros(1,length(theta_start:d_theta:theta_end));

counter = 1;

for theta = theta_start:d_theta:theta_end

G_t2(counter) = func(-theta);

counter = counter + 1;

end

theta = theta_start:d_theta:theta_end;


L_p_1 = (4*pi*lambda/R)^2;


figure;

plot(theta,20*log10(G_t.*G_t2*L_p_1/F));

xlim([theta_start theta_end])

function F = func(theta)

  fun = cos((pi/2)*cos(theta))/sin(theta);

  F = abs(fun)^2;

end



