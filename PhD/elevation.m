clc;
clear;
close all;


%% data load

H_0 = 1288;

h_0 = 815; %m

avg_terrain_height = H_0-h_0;

r = 6371e3; %m

reference = 1e3; %m

A = readtable('Blacksburg_elevation.txt');

A2 = table2array(A(:,2))*pi/180; %lat, rad

A3 = table2array(A(:,3))*pi/180; %long, rad

A4 = table2array(A(:,4)); %elevation, m

stable_A4 = A4; % elevation above sea level, m

f1 = 1256e6; %MHz

f2 = 1292e6; %MHz

c = 3e8; %m/s

lambda1 = c/f1;

lambda2 = c/f2;

R_b1 = zeros(1,length(A3));

R_b2 = zeros(1,length(A3));

M = median(A4);

for counter = 1:length(A3)

    %h_1 = HAAT(A2(counter), A3(counter), A4(counter), r, reference);

    h_1 = (A4(counter)-min(A4))+1;

    R_b1(counter) = abs(4*pi/lambda1)*(H_0-min(A4))*h_1;

    R_b2(counter) = abs(4*pi/lambda2)*(H_0-min(A4))*h_1;

end

points = 500;

phi_1 = 37.517*pi/180; %rad base station

theta_1 = -79.510*pi/180; %rad


phi_pem = 37.3196; %deg pembroke

theta_pem = -80.6390; %deg


phi_roa = 37.2710; %deg roanoke

theta_roa = -79.9414; %deg


theta_b = -80.4139; %deg, bb
    
phi_b = 37.2296;% deg


theta_c = -80.4089; %deg, c
    
phi_c = 37.1299;% deg

x_lower = -80.75;

x_upper = -78.75;

y_lower = 36.75;

y_upper = 38.2;

clim_lower = 0;

clim_upper = 200;

R = zeros(1,length(A2));

for counter = 1:length(A3)

    Delta_phi = A2(counter)-phi_1;

    Delta_theta = A3(counter)-theta_1;

    a = sin(Delta_phi/2)^2+cos(phi_1)*cos(A2(counter))*sin(Delta_theta/2)^2;

    c = 2*atan2(sqrt(a),sqrt(1-a));

    R(counter) = r*c;

end


%% for n = 2

n = 2;

L_p1_n2 = zeros(1,length(A3));

L_p2_n2 = zeros(1,length(A3));

for counter = 1:length(A3)

    L_p1_n2(counter) = path_loss(lambda1, n, R(counter), R_b1(counter), h_0, A4(counter));

    L_p2_n2(counter) = path_loss(lambda2, n, R(counter), R_b2(counter), h_0, A4(counter));

end


%% plotting path loss for n = 2

% Define grid range
yRange = linspace(min(A2), max(A2), points); % 100 points for better resolution
xRange = linspace(min(A3), max(A3), points); % 100 points for better resolution

% Create grid
[XGrid, YGrid] = meshgrid(xRange, yRange);

plotting(A3, A2, L_p1_n2, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f1)


plotting(A3, A2, L_p2_n2, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f2)




%% for n = 5

n = 5;

L_p1_n5 = zeros(1,length(A3));

L_p2_n5 = zeros(1,length(A3));

for counter = 1:length(A3)

    L_p1_n5(counter) = path_loss(lambda1, n, R(counter), R_b1(counter), h_0, A4(counter));

    L_p2_n5(counter) = path_loss(lambda2, n, R(counter), R_b2(counter), h_0, A4(counter));

end

%% 

plotting(A3, A2, L_p1_n5, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f1)


hold all;



plotting(A3, A2, L_p2_n5, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f2)

hold all

%% heat map for elevation

ZGrid = griddata(A3, A2, stable_A4, XGrid, YGrid, 'natural');

% Display the grid with imagesc
figure;
imagesc(xRange*180/pi, yRange*180/pi, ZGrid);

colormap("jet")
grid on;
ylabel('Latitude [deg]')
xlabel('Longitude [deg]')
set (gca,'YDir','normal')
set(gca,'FontSize',20)

axis('square')

hc = colorbar('EastOutside');

ylabel(hc,'Elevation [m]','FontSize',16)

title('Elevation Map')

hold all;

axis('square')

hc = colorbar('EastOutside');

ylabel(hc,'Elevation [m]','FontSize',16)

hold all;

hold on

plot(theta_1*180/pi,phi_1*180/pi, 'MarkerFaceColor', [1 0 0], 'Marker', '^', 'MarkerSize', 15);

hold all;

plot(theta_b, phi_b,'MarkerFaceColor', [0 0 0], 'Marker', '^', 'MarkerSize', 15);

plot(theta_c,phi_c, 'MarkerFaceColor', [0.5 0 0.3], 'Marker', '^', 'MarkerSize', 15);

plot(theta_pem,phi_pem, 'MarkerFaceColor', [0.3 0 0.5], 'Marker', '^', 'MarkerSize', 15);

plot(theta_roa,phi_roa, 'MarkerFaceColor', [0.5 0.5 0.3], 'Marker', '^', 'MarkerSize', 15);

legend('Bedford, VA (Base Station)','Blacksburg, VA', 'Christiansburg, VA', 'Pembroke, VA', 'Roanoke, VA')

xlim([x_lower x_upper])
ylim([y_lower y_upper])

hold all;

%% test

tx = txsite(Latitude= phi_1*180/pi,Longitude= theta_1*180/pi, ...
    TransmitterFrequency=f1);
rx = rxsite(Latitude=phi_b,Longitude= theta_b);


pm = propagationModel("freespace");
pl = pathloss(pm,rx,tx);
disp(pl)


%% functions

function L = path_loss(lambda, n, R, R_b, h , h_0)

R_0 = 4.12*sqrt(h)*1000; %radio horizon distance

%if(R<=R_0 || h>=h_0)
%if(1==1)
%if(h >= h_0)
if(R<=R_0)
%if(R<=R_0 && h>=h_0)

    


    if ((R)>(R_b))

        L_p = (lambda/(4*pi*R_b))^(-2);

        L =  10*log10(abs(L_p...
            *(R/R_b)^n));

    else

        L_p = (lambda/(4*pi*R))^(-2);

        L = 10*log10(abs(L_p));

    end

else

    L = Inf; %dB

end

end


function plotting(A3, A2, L_p2_n5, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f)

    figure;

    hold all;

    ZGrid = griddata(A3, A2, L_p2_n5, XGrid, YGrid, 'natural');
    
    % Display the grid with imagesc
    imagesc(xRange*180/pi, yRange*180/pi, ZGrid);
    
    colormap("jet")
    grid on;
    ylabel('Latitude [deg]')
    xlabel('Longitude [deg]')
    set (gca,'YDir','normal')
    set(gca,'FontSize',20)
    
    axis('square')
    
    hc = colorbar('EastOutside');
    
    ylabel(hc,'Pathloss [dB]','FontSize',16)

    Ab = ['f = ', num2str(f/(1e6)), ' MHz', ', n = ', num2str(n)];
    
    title(Ab)
    
    hold all;
    
    hold on
    
    plot(theta_1*180/pi,phi_1*180/pi, 'MarkerFaceColor', [1 0 0], 'Marker', '^', 'MarkerSize', 15);
    
    hold all;
    
    plot(theta_b, phi_b,'MarkerFaceColor', [0 0 0], 'Marker', '^', 'MarkerSize', 15);
    
    plot(theta_c, phi_c, 'MarkerFaceColor', [0.5 0 0.3], 'Marker', '^', 'MarkerSize', 15);
    
    plot(theta_pem,phi_pem, 'MarkerFaceColor', [0.3 0 0.5], 'Marker', '^', 'MarkerSize', 15);
    
    plot(theta_roa,phi_roa, 'MarkerFaceColor', [0.5 0.5 0.3], 'Marker', '^', 'MarkerSize', 15);
    
    legend('Bedford, VA (Base Station)','Blacksburg, VA', 'Christiansburg, VA', 'Pembroke, VA', 'Roanoke, VA')
    
    xlim([x_lower x_upper])
    ylim([y_lower y_upper])
    
    clim([clim_lower clim_upper])

end

function A4 = HAAT(A2, A3, A4, r, reference)

temp_h = zeros(1,length(A2));

for counter = 1:length(A3)

    ref_lat = A2(counter);

    ref_long = A3(counter);

    sum = 0;

    index = 0;

    for counter2 = 1:length(A3)

        Delta_phi = A2(counter2)-ref_lat;

        Delta_theta = A3(counter2)-ref_long;

        a = sin(Delta_phi/2)^2+cos(ref_lat)*cos(A2(counter))*sin(Delta_theta/2)^2;

        c = 2*atan2(sqrt(a),sqrt(1-a));

        dist = r*c;

        if (dist <= reference)

            sum = sum + A4(counter2);

            index = index + 1;

        end

        temp_h(counter) = sum/index;%

    end

end

A4 = A4-temp_h';

end
