clc;
clear;
close all;

f = 1e9;

c = 3e8;

lambda = c/f;

beta = 2*pi/lambda;

a = 60*lambda;

omega = 2*pi*f;

epsilon = 8.85e-12;

mu = (4*pi)*(10^-7);

eta = sqrt(mu/epsilon);

%% outward power radiated by dipole

sum = 0;

theta = 0;

d_theta = pi/90;

d_x = lambda/10;

d_y = d_x;

d_aperture = d_x*d_y;

d_phi = d_theta;

while (theta <= pi)

    phi = 0;

    while (phi<=2*pi)

        fun = (-(cos(phi)^2)*(sin(theta)^3)+sin(theta))*d_theta*d_phi;
        sum = sum + fun;
        phi = phi + d_phi;

    end

    theta = theta + d_theta;
end


term = ((omega*mu)^2)/(128*eta*pi^2);

P_dipole = term*sum;



%% power radiated by the aperture

%calculating M_s

z = 100*lambda;

x_length = length(0:d_x:a);

y_length = length(0:d_y:a);

counter = 1;

for x = 0:d_x:a

    for y=0:d_y:a

        if (sqrt(x^2+y^2) <= a)

            counter = counter + 1;

        end

    end
end

M_s  = zeros(counter,3);

r_hat_prime = zeros(counter,3);

n_hat = [0 0 1];

counter = 1;

for x = 0:d_x:a

    for y=0:d_y:a

        if (sqrt(x^2+y^2) <= a)

            [r,theta,phi] = spherical(x,y,z);

            theta_hat = theta_h(theta,phi);

            phi_hat = phi_h(theta,phi);

            polarization = theta_hat*cos(theta)*cos(phi)-phi_hat*sin(phi);

            greens = exp(-1i*beta*r)/r;

            constants = -1i*omega*mu/(8*pi);

            E_a = constants*greens*polarization;

            M_s(counter,:) = -2*cross(n_hat,E_a);

            r_hat_prime(counter,:) = [x y 0];

            counter = counter + 1;


        end



    end


end


%% check

P_c = 0;

for x = 0:d_x:a

    for y=0:d_y:a

        if (sqrt(x^2+y^2) <= a)

            [r,theta,phi] = spherical(x,y,z);

            theta_hat = theta_h(theta,phi);

            phi_hat = phi_h(theta,phi);

            polarization = theta_hat*cos(theta)*cos(phi)-phi_hat*sin(phi);

            greens = exp(-1i*beta*r)/r;

            E_a = constants*greens*polarization;

            P_c = P_c + (norm(E_a)^2)*d_aperture;


        end



    end


end

P_c = P_c/(2*eta);

%%

%calculating F

x = 0:d_x:a;

y = 0:d_y:a;

counter = 1;

sum = 0;

for theta = 0:d_theta:pi/2

    for phi=0:d_phi:2*pi

        r_hat = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

        theta_hat = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];

        phi_hat = [-sin(phi) cos(phi) 0];

        H = 0;

        counter = 1;

        for x = 0:d_x:a

            for y=0:d_y:a

                if (sqrt(x^2+y^2) <= a)

                    product = dot(r_hat,r_hat_prime(counter,:));

                    integrand = (epsilon/(4*pi))*M_s(counter,:)*exp(1i*beta*product)*d_aperture;

                    H = H + integrand;

                    counter = counter + 1;

                end
            end
        end

        H = -1i*omega*(theta_hat*dot(theta_hat,H) + phi_hat*dot(phi_hat,H));

        counter = 1;

        E = eta*cross(H,r_hat);

%         E = 0;
%
%         for x = 0:d_x:a
% 
%             for y=0:d_y:a
% 
%                 if (sqrt(x^2+y^2) <= a)
% 
%                     phasing = phi_hat*dot(theta_hat,M_s(counter,:)) -theta_hat*dot(phi_hat,M_s(counter,:));
% 
%                     phasing = cross(phasing, r_hat);
% 
%                     product = dot(r_hat,r_hat_prime(counter,:));
% 
%                     integrand = (epsilon*eta/(4*pi))*phasing*exp(1i*beta*product)*d_aperture;
% 
%                     E = E + integrand;
% 
%                     counter = counter + 1;
% 
%                 end
% 
%             end
% 
% 
%         end

        integrand = cross(E,conj(H));

        sum = sum+dot(r_hat,integrand*sin(theta)*d_phi*d_theta);

        %sum = sum+(1/(2*eta))*(norm(E)^2)*sin(theta)*d_phi*d_theta;

    end


end

P_radiated = (1/2)*(real(sum));










