function [a_p, a_m] = solve_a(n, phi, phi_p, pos_min, N_p, N_m)

if(pos_min == 1)

    phase = phi + phi_p;

else
    phase = phi - phi_p;

end

num_p = 2*n*pi*N_p - phase;

den_p = 2;

a_p = 2*cos(num_p/den_p)^2;


num_m = 2*n*pi*N_m - phase;

den_m = 2;

a_m = 2*cos(num_m/den_m)^2;



end