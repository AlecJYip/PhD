clc;
clear;
close all;


%% data load

A = readtable('blacksburg.txt');

A2 = table2array(A(:,2))*pi/180; %lat, rad

A3 = table2array(A(:,3))*pi/180; %long, rad

A4 = table2array(A(:,4));%*12*2.54/100; %elevation, m

H_0 = 1288;

% h_0 = 815; %m

%h_0 = H_0-min(A4);

r = 6371e3; %m

reference = 20e3; %m

h_0 = H_0-sum(A4)/length(A4);

avg_terrain_height = H_0-h_0;

stable_A4 = A4; % elevation above sea level, m

f1 = 1215e6; %MHz

f2 = 1350e6; %MHz

c = 3e8; %m/s

lambda1 = c/f1;

lambda2 = c/f2;

R_b1 = zeros(1,length(A3));

R_b2 = zeros(1,length(A3));

M = median(A4);



points = 500;

phi_1 = 37.517*pi/180; %rad base station

theta_1 = -79.510*pi/180; %rad


phi_pem = 37.3196; %deg pembroke

theta_pem = -80.6390; %deg

elev_pem = 502.67; %m


%(37.3196, -80.6390)


phi_roa = 37.2710; %deg roanoke

theta_roa = -79.9414; %deg

elev_roa = 284.24; %m


theta_b = -80.4139; %deg, bb
    
phi_b = 37.2296;% deg

elev_bb = 632.00; %m


theta_c = -80.4089; %deg, c
    
phi_c = 37.1299;% deg

elev_c = 635.01; %m



x_lower = -80.75;

x_upper = -78.75;

y_lower = 36.75;

y_upper = 38.2;


for counter = 1:length(A3)

    R_b1(counter) = ((4*pi/lambda1)*(h_0)*(A4(counter)-avg_terrain_height));

    R_b2(counter) = ((4*pi/lambda2)*(h_0)*(A4(counter)-avg_terrain_height));

    if(R_b1(counter) < 0)

        R = dist(A2(counter), A3(counter), phi_1,theta_1,r);

        R_b1(counter) = A4(counter)*R/(H_0+A4(counter));

        R_b2(counter) = R_b1(counter);

    else

        R = dist(A2(counter), A3(counter), phi_1,theta_1,r);

        R_b1(counter) = h_0*R/(H_0+A4(counter));

        R_b2(counter) = R_b1(counter);

    end

end

% x_lower = -82;
% 
% x_upper = -78;
% 
% y_lower = 36.5;
% 
% y_upper = 38;

clim_lower = 0;

clim_upper = 200;

clear R;

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

figure;

subplot(1,2,1)

plotting(A3, A2, L_p1_n2, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f1)

hold all

subplot(1,2,2)

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


figure;

subplot(1, 2, 1)

plotting(A3, A2, L_p1_n5, XGrid, YGrid, theta_1, phi_1,...
    theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
    x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
    xRange, yRange, n, f1)


hold all;

%legend off

%colorbar off

subplot(1,2,2)

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

%% Individual Path Loss (that we calculate) for n = 2

n = 2;


disp('For 1256 MHz')

R = dist(phi_b*pi/180, theta_b*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_b*pi/180, theta_b*pi/180, elev_bb,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_bb);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Blacksburg: ', num2str(pl_n_bb), ' dB'];
disp(a)


R = dist(phi_c*pi/180, theta_c*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_c*pi/180, theta_c*pi/180, elev_c,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_c);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Christiansburg: ', num2str(pl_n_bb), ' dB'];
disp(a)



R = dist(phi_roa*pi/180, theta_roa*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_roa*pi/180, theta_roa*pi/180, elev_roa,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_roa);


a = ['Calculated Pathloss for n = ', num2str(n), ' at Roanoke: ', num2str(pl_n_bb), ' dB'];
disp(a)

R = dist(phi_pem*pi/180, theta_pem*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_pem*pi/180, theta_pem*pi/180, elev_pem,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_pem);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Pembroke: ', num2str(pl_n_bb), ' dB'];
disp(a)


L_p_2_power_f2 =[];

disp('For 1292 MHz')

R = dist(phi_b*pi/180, theta_b*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_b*pi/180, theta_b*pi/180, elev_bb,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_bb);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Blacksburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_2_power_f2 = [L_p_2_power_f2 pl_n_bb];


R = dist(phi_c*pi/180, theta_c*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_c*pi/180, theta_c*pi/180, elev_c,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_c);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Christiansburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_2_power_f2 = [L_p_2_power_f2 pl_n_bb];



R = dist(phi_roa*pi/180, theta_roa*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_roa*pi/180, theta_roa*pi/180, elev_roa,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_roa);


a = ['Calculated Pathloss for n = ', num2str(n), ' at Roanoke: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_2_power_f2 = [L_p_2_power_f2 pl_n_bb];



R = dist(phi_pem*pi/180, theta_pem*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_pem*pi/180, theta_pem*pi/180, elev_pem,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_pem);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Pembroke: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_2_power_f2 = [L_p_2_power_f2 pl_n_bb];


%% Individual Path Loss (that we calculate) for n = 5

n = 5;

L_p_5_power_f1 =[];


disp('For 1256 MHz')

R = dist(phi_b*pi/180, theta_b*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_b*pi/180, theta_b*pi/180, elev_bb,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_bb);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Blacksburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f1 = [L_p_5_power_f1 pl_n_bb];



R = dist(phi_c*pi/180, theta_c*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_c*pi/180, theta_c*pi/180, elev_c,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_c);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Christiansburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f1 = [L_p_5_power_f1 pl_n_bb];






R = dist(phi_roa*pi/180, theta_roa*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_roa*pi/180, theta_roa*pi/180, elev_roa,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_roa);
L_p_5_power_f1 = [L_p_5_power_f1 pl_n_bb];



a = ['Calculated Pathloss for n = ', num2str(n), ' at Roanoke: ', num2str(pl_n_bb), ' dB'];
disp(a)


R = dist(phi_pem*pi/180, theta_pem*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda1, h_0, H_0, phi_pem*pi/180, theta_pem*pi/180, elev_pem,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda1, n, R, R_b1, h_0, elev_pem);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Pembroke: ', num2str(pl_n_bb), ' dB'];
disp(a)
L_p_5_power_f1 = [L_p_5_power_f1 pl_n_bb];


disp('For 1292 MHz')

L_p_5_power_f2 = [];


R = dist(phi_b*pi/180, theta_b*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_b*pi/180, theta_b*pi/180, elev_bb,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_bb);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Blacksburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f2 = [L_p_5_power_f2 pl_n_bb];



R = dist(phi_c*pi/180, theta_c*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_c*pi/180, theta_c*pi/180, elev_c,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_c);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Christiansburg: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f2 = [L_p_5_power_f2 pl_n_bb];




R = dist(phi_roa*pi/180, theta_roa*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_roa*pi/180, theta_roa*pi/180, elev_roa,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_roa);


a = ['Calculated Pathloss for n = ', num2str(n), ' at Roanoke: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f2 = [L_p_5_power_f2 pl_n_bb];


R = dist(phi_pem*pi/180, theta_pem*pi/180, phi_1,theta_1,r);
Aa = break_point(lambda2, h_0, H_0, phi_pem*pi/180, theta_pem*pi/180, elev_pem,...
    phi_1,theta_1, avg_terrain_height, r);
R_b1 = Aa(1);
pl_n_bb = path_loss(lambda2, n, R, R_b1, h_0, elev_pem);
a = ['Calculated Pathloss for n = ', num2str(n), ' at Pembroke: ', num2str(pl_n_bb), ' dB'];
disp(a)

L_p_5_power_f2 = [L_p_5_power_f2 pl_n_bb];

%% test



tx = txsite(Latitude= phi_1*180/pi,Longitude= theta_1*180/pi, ...
    TransmitterFrequency=f1);
pm = propagationModel("freespace");

L_p_2_power_f1 = [];

pl = display_pathloss(pm, theta_b, phi_b, tx);

L_p_2_power_f1 = [L_p_2_power_f1 pl];

a  = ['Freespace Pathloss (1256 MHz) for Blacksburg, VA: ', num2str(pl),' dB'];
disp(a)


pl = display_pathloss(pm, theta_c, phi_c, tx);

L_p_2_power_f1 = [L_p_2_power_f1 pl];


a  = ['Freespace Pathloss (1256 MHz) for Christiansburg, VA: ', num2str(pl),' dB'];
disp(a)


pl = display_pathloss(pm, theta_roa, phi_roa, tx);

L_p_2_power_f1 = [L_p_2_power_f1 pl];


a  = ['Freespace Pathloss (1256 MHz) for Roanoke, VA: ', num2str(pl),' dB'];
disp(a)


pl = display_pathloss(pm, theta_pem, phi_pem, tx);

L_p_2_power_f1 = [L_p_2_power_f1 pl];


a  = ['Freespace Pathloss (1256 MHz) for Pembroke, VA: ', num2str(pl),' dB'];
disp(a)



pm = propagationModel("longley-rice");

pl = display_pathloss(pm, theta_roa, phi_roa, tx);

L_p_longley_rice_f1 = [];

L_p_longley_rice_f1 = [L_p_longley_rice_f1 pl];

a  = ['Longley-Rice Pathloss (1256 MHz) for Blacksburg, VA: ', num2str(pl),' dB'];
disp(a)

pl = display_pathloss(pm, theta_c, phi_c, tx);

a  = ['Longley-Rice Pathloss (1256 MHz) for Christiansburg, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f1 = [L_p_longley_rice_f1 pl];





pl = display_pathloss(pm, theta_roa, phi_roa, tx);

a  = ['Longley-Rice Pathloss (1256 MHz) for Roanoke, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f1 = [L_p_longley_rice_f1 pl];

pl = display_pathloss(pm, theta_pem, phi_pem, tx);

a  = ['Longley-Rice Pathloss (1256 MHz) for Pembroke, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f1 = [L_p_longley_rice_f1 pl];



tx = txsite(Latitude= phi_1*180/pi,Longitude= theta_1*180/pi, ...
    TransmitterFrequency=f2);
pm = propagationModel("freespace");


pl = display_pathloss(pm, theta_b, phi_b, tx);

a  = ['Freespace Pathloss (1292 MHz) for Blacksburg, VA: ', num2str(pl),' dB'];
disp(a)


pl = display_pathloss(pm, theta_c, phi_c, tx);

a  = ['Freespace Pathloss (1292 MHz) for Christiansburg, VA: ', num2str(pl),' dB'];
disp(a)



pl = display_pathloss(pm, theta_roa, phi_roa, tx);

a  = ['Freespace Pathloss (1292 MHz) for Roanoke, VA: ', num2str(pl),' dB'];
disp(a)



pl = display_pathloss(pm, theta_pem, phi_pem, tx);

a  = ['Freespace Pathloss (1292 MHz) for Pembroke, VA: ', num2str(pl),' dB'];
disp(a)




L_p_longley_rice_f2 = [];

pm = propagationModel("longley-rice");

pl = display_pathloss(pm, theta_b, phi_b, tx);


a  = ['Longley-Rice Pathloss (1292 MHz) for Blacksburg, VA: ', num2str(pl),' dB'];
disp(a)

pl = display_pathloss(pm, theta_c, phi_c, tx);

L_p_longley_rice_f2 = [L_p_longley_rice_f2 pl];


a  = ['Longley-Rice Pathloss (1292 MHz) for Christiansburg, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f2 = [L_p_longley_rice_f2 pl];



pl = display_pathloss(pm, theta_roa, phi_roa, tx);

a  = ['Longley-Rice Pathloss (1292 MHz) for Roanoke, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f2 = [L_p_longley_rice_f2 pl];


pl = display_pathloss(pm, theta_pem, phi_pem, tx);

a  = ['Longley-Rice Pathloss (1292 MHz) for Pembroke, VA: ', num2str(pl),' dB'];
disp(a)

L_p_longley_rice_f2 = [L_p_longley_rice_f2 pl];



%% Longley-Rice Coverage
% 
% 
% tx = txsite(Latitude= phi_1*180/pi,Longitude= theta_1*180/pi, ...
%     TransmitterFrequency=f1);
% pm = propagationModel("longley-rice");
% 
% L_longley_p1 = zeros(1,length(A3));
% 
% for counter = 1:length(A3)
% 
%     L_longley_p1(counter) = display_pathloss(pm, A3(counter)*180/pi, A2(counter)*180/pi, tx);
%     %disp(counter)
% 
% end
% 
% figure;
% 
% subplot(1,2,1);
% 
% hold all;
% 
% 
% plotting(A3, A2, L_longley_p1, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz'];
% title(a);
% hold all
% 
% subplot(1,2,2)
% 
% hold all
% 
% tx = txsite(Latitude= phi_1*180/pi,Longitude= theta_1*180/pi, ...
%     TransmitterFrequency=f2);
% pm = propagationModel("longley-rice");
% 
% L_longley_p2 = zeros(1,length(A3));
% 
% for counter = 1:length(A3)
% 
%     L_longley_p2(counter) = display_pathloss(pm, A3(counter)*180/pi, A2(counter)*180/pi, tx);
%     %disp(counter)
% 
% end
% plotting(A3, A2, L_longley_p2, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% a = ['f = ', num2str(f2/(1e6)), ' MHz'];
% title(a);
% hold all


%% Big plot of PL for f = f1


% figure;
% 
% subplot(2,2,1);
% 
% hold all;
% 
% n = 2;
% 
% plotting(A3, A2, L_p1_n2, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law for n = ', num2str(n)];
% title(a);
% hold all
% 
% legend off
% 
% subplot(2,2,3);
% 
% hold all;
% 
% n = 5;
% 
% plotting(A3, A2, L_p1_n5, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law for n = ', num2str(n)];
% title(a);
% hold all
% 
% legend off
% 
% subplot(2,2,4);
% 
% hold all;
% 
% 
% plotting(A3, A2, L_longley_p1, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all


%% Big plot of PL for f = f2


% figure;
% 
% subplot(2,2,1);
% 
% hold all;
% 
% n = 2;
% 
% plotting(A3, A2, L_p2_n2, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% title off
% hold all
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law for n = ', num2str(n)];
% title(a);
% hold all
% 
% legend off
% 
% subplot(2,2,3);
% 
% hold all;
% 
% n = 5;
% 
% plotting(A3, A2, L_p2_n5, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% title off
% hold all
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law for n = ', num2str(n)];
% title(a);
% hold all
% 
% legend off
% 
% subplot(2,2,4);
% 
% hold all;
% 
% 
% plotting(A3, A2, L_longley_p2, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% title off
% hold all
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all

%% Input SNR

G_T = 10^(35/10); %linear gain of transmitting antenna

P_T = 2e6; %watts, transmit power

F = 9; %dB

G_R = 10^(8/10); %linear gain of receiver antenna

B_RX = 2.4e6; %bandwidth of receiver

T_RX = 290; %kelvin

k = 1.38e-23; %J/s, boltzmanns constant

a = load('temp_data.txt');

epsilon_R = 1;

eta_R = 1;

epsilon_T = 1;

eta_T = 1;

epsilon_p = 0.5;

bb_high = a(1,:)+273.15; %kelvin

bb_low = a(2,:)+273.15; %kelvin

c_high = a(3,:)+273.15; %kelvin

c_low = a(4,:)+273.15; %kelvin

roa_high = a(5,:)+273.15; %kelvin

roa_low = a(6,:)+273.15; %kelvin

pem_high = a(7,:)+273.15; %kelvin

pem_low = a(8,:)+273.15; %kelvin


%% for n = 2



SNR_high_bb_n_2 = zeros(1,length(bb_low));

SNR_high_c_n_2 = zeros(1,length(bb_low));

SNR_high_roa_n_2 = zeros(1,length(bb_low));

SNR_high_pem_n_2 = zeros(1,length(bb_low));

SNR_low_bb_n_2 = zeros(1,length(bb_low));

SNR_low_c_n_2 = zeros(1,length(bb_low));

SNR_low_roa_n_2 = zeros(1,length(bb_low));

SNR_low_pem_n_2 = zeros(1,length(bb_low));


for counter = 1:length(L_p_longley_rice_f1)

    L_p_longley_rice_f1(counter) = 10^(L_p_longley_rice_f1(counter)/10);

    L_p_2_power_f1(counter) = 10^(L_p_2_power_f1(counter)/10);

    L_p_5_power_f1(counter) = 10^(L_p_5_power_f1(counter)/10);

    L_p_longley_rice_f2(counter) = 10^(L_p_longley_rice_f2(counter)/10);

    L_p_2_power_f2(counter) = 10^(L_p_2_power_f2(counter)/10);

    L_p_5_power_f2(counter) = 10^(L_p_5_power_f2(counter)/10);

end


for counter = 1:length(SNR_low_pem_n_2)

    SNR_high_bb_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f1(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

end

figure; 

subplot(2,2,1);


hold all;

plot(SNR_high_bb_n_2, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_n_2, 'r--', 'LineWidth',2); 

plot(SNR_high_c_n_2, 'b-', 'LineWidth',2); 

plot(SNR_low_c_n_2, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_n_2, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_n_2, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_n_2, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_n_2, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f1/(1e6)), ' MHz, Power Law Model for n = 2'];

title(a)

%legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
%     'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
%     'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
%     'Pembroke (Avg High)', 'Pembroke (Avg Low)')


%for n = 5

SNR_high_bb_n_5 = zeros(1,length(bb_low));

SNR_high_c_n_5 = zeros(1,length(bb_low));

SNR_high_roa_n_5 = zeros(1,length(bb_low));

SNR_high_pem_n_5 = zeros(1,length(bb_low));

SNR_low_bb_n_5 = zeros(1,length(bb_low));

SNR_low_c_n_5 = zeros(1,length(bb_low));

SNR_low_roa_n_5 = zeros(1,length(bb_low));

SNR_low_pem_n_5 = zeros(1,length(bb_low));


for counter = 1:length(SNR_low_pem_n_5)

    SNR_high_bb_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f1(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

end
xlim([1 12])

subplot(2,2,3); hold all;

plot(SNR_high_bb_n_5, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_n_5, 'r--', 'LineWidth',2); 

plot(SNR_high_c_n_5, 'b-', 'LineWidth',2); 

plot(SNR_low_c_n_5, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_n_5, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_n_5, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_n_5, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_n_5, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f1/(1e6)), ' MHz, Power Law Model for n=5'];

title(a)

% legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
%     'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
%     'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
%     'Pembroke (Avg High)', 'Pembroke (Avg Low)')

xlim([1 12])
%Longley-Rice Model


SNR_high_bb_lr = zeros(1,length(bb_low));

SNR_high_c_lr = zeros(1,length(bb_low));

SNR_high_roa_lr = zeros(1,length(bb_low));

SNR_high_pem_lr = zeros(1,length(bb_low));

SNR_low_bb_lr = zeros(1,length(bb_low));

SNR_low_c_lr = zeros(1,length(bb_low));

SNR_low_roa_lr = zeros(1,length(bb_low));

SNR_low_pem_lr = zeros(1,length(bb_low));


for counter = 1:length(SNR_low_pem_lr)

    SNR_high_bb_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f1(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);


end


subplot(2,2,4) ; hold all;

plot(SNR_high_bb_lr, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_lr, 'r--', 'LineWidth',2); 

plot(SNR_high_c_lr, 'b-', 'LineWidth',2); 

plot(SNR_low_c_lr, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_lr, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_lr, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_lr, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_lr, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f1/(1e6)), ' MHz, Longley-Rice Model'];

title(a)

legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
    'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
    'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
    'Pembroke (Avg High)', 'Pembroke (Avg Low)')
xlim([1 12])


%% for n = 5



SNR_high_bb_n_2 = zeros(1,length(bb_low));

SNR_high_c_n_2 = zeros(1,length(bb_low));

SNR_high_roa_n_2 = zeros(1,length(bb_low));

SNR_high_pem_n_2 = zeros(1,length(bb_low));

SNR_low_bb_n_2 = zeros(1,length(bb_low));

SNR_low_c_n_2 = zeros(1,length(bb_low));

SNR_low_roa_n_2 = zeros(1,length(bb_low));

SNR_low_pem_n_2 = zeros(1,length(bb_low));


for counter = 1:length(SNR_low_pem_n_2)

    SNR_high_bb_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_n_2(counter) = input_SNR(P_T, G_R, G_T, L_p_2_power_f2(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

end

figure;

subplot(2,2,1);

hold all;

plot(SNR_high_bb_n_2, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_n_2, 'r--', 'LineWidth',2); 

plot(SNR_high_c_n_2, 'b-', 'LineWidth',2); 

plot(SNR_low_c_n_2, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_n_2, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_n_2, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_n_2, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_n_2, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f2/(1e6)), ' MHz, Power Law Model for n = 2'];

title(a)

% legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
%     'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
%     'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
%     'Pembroke (Avg High)', 'Pembroke (Avg Low)')


%for n = 5

SNR_high_bb_n_5 = zeros(1,length(bb_low));

SNR_high_c_n_5 = zeros(1,length(bb_low));

SNR_high_roa_n_5 = zeros(1,length(bb_low));

SNR_high_pem_n_5 = zeros(1,length(bb_low));

SNR_low_bb_n_5 = zeros(1,length(bb_low));

SNR_low_c_n_5 = zeros(1,length(bb_low));

SNR_low_roa_n_5 = zeros(1,length(bb_low));

SNR_low_pem_n_5 = zeros(1,length(bb_low));


for counter = 1:length(SNR_low_pem_n_5)

    SNR_high_bb_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_n_5(counter) = input_SNR(P_T, G_R, G_T, L_p_5_power_f2(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

end
xlim([1 12])

subplot(2,2,3); hold all;

plot(SNR_high_bb_n_5, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_n_5, 'r--', 'LineWidth',2); 

plot(SNR_high_c_n_5, 'b-', 'LineWidth',2); 

plot(SNR_low_c_n_5, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_n_5, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_n_5, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_n_5, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_n_5, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f2/(1e6)), ' MHz, Power Law Model for n=5'];

title(a)

% legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
%     'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
%     'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
%     'Pembroke (Avg High)', 'Pembroke (Avg Low)')

xlim([1 12])
%Longley-Rice Model


SNR_high_bb_lr = zeros(1,length(bb_low));

SNR_high_c_lr = zeros(1,length(bb_low));

SNR_high_roa_lr = zeros(1,length(bb_low));

SNR_high_pem_lr = zeros(1,length(bb_low));

SNR_low_bb_lr = zeros(1,length(bb_low));

SNR_low_c_lr = zeros(1,length(bb_low));

SNR_low_roa_lr = zeros(1,length(bb_low));

SNR_low_pem_lr = zeros(1,length(bb_low));


for counter = 1:length(SNR_low_pem_lr)

    SNR_high_bb_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(1), k, bb_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_c_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(2), k, c_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_roa_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(3), k, roa_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_high_pem_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(4), k, pem_high(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_bb_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(1), k, bb_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_c_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(2), k, c_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_roa_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(3), k, roa_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);

    SNR_low_pem_lr(counter) = input_SNR(P_T, G_R, G_T, L_p_longley_rice_f2(4), k, pem_low(counter), T_RX, B_RX, epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F);


end


subplot(2,2,4); hold all;

plot(SNR_high_bb_lr, 'r-', 'LineWidth',2); 

plot(SNR_low_bb_lr, 'r--', 'LineWidth',2); 

plot(SNR_high_c_lr, 'b-', 'LineWidth',2); 

plot(SNR_low_c_lr, 'b--', 'LineWidth',2); 

plot(SNR_high_roa_lr, 'color', '#f3b994', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_roa_lr, 'color', '#f3b994', 'LineStyle','--','LineWidth',2); 

plot(SNR_high_pem_lr, 'color', '#90EE90', 'LineStyle','-','LineWidth',2); 

plot(SNR_low_pem_lr, 'color', '#90EE90', 'LineStyle','--','LineWidth',2); 

grid on;

ylabel('SNR [dB]')

xlabel('Month')

set(gca,'FontSize', 20)

axis square

a = ['f = ',num2str(f2/(1e6)), ' MHz, Longley-Rice Model'];

title(a)

legend('Blacksburg (Avg High)', 'Blacksburg (Avg Low)',...
    'Christiansburg (Avg High)', 'Christiansburg, (Avg Low)',...
    'Roanoke (Avg High)', 'Roanoke (Avg Low)',...
    'Pembroke (Avg High)', 'Pembroke (Avg Low)')
xlim([1 12])


%% received power 

% for counter = 1:length(L_p1_n5)
% 
%     L_p1_n2(counter) = 10^(0.1*L_p1_n2(counter));
% 
%     L_p2_n2(counter) = 10^(0.1*L_p2_n2(counter));
% 
%     L_p1_n5(counter) = 10^(0.1*L_p1_n5(counter));
% 
%     L_p2_n5(counter) = 10^(0.1*L_p2_n5(counter));
% 
%     L_longley_p1(counter) = 10^(0.1*L_longley_p1(counter));
% 
%     L_longley_p2(counter) = 10^(0.1*L_longley_p2(counter));
% 
% 
% end

% %%
% 
% 
% P_R_p1_n2 = received_power(P_T, G_T, G_R, L_p1_n2);
% 
% P_R_p2_n2 = received_power(P_T, G_T, G_R, L_p2_n2);
% 
% P_R_p1_n5 = received_power(P_T, G_T, G_R, L_p1_n5);
% 
% P_R_p2_n5 = received_power(P_T, G_T, G_R, L_p2_n5);
% 
% P_longley_p1 = received_power(P_T, G_T, G_R, L_longley_p1);
% 
% P_longley_p2 = received_power(P_T, G_T, G_R, L_longley_p2);
% 
% 
% 
% %% plotting power
% 
% 
% figure;
% 
% subplot(2,2,1);
% 
% hold all;
% 
% n = 2;
% 
% plotting(A3, A2, P_R_p1_n2+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% 
% legend off
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% legend off
% 
% 
% subplot(2,2,3)
% 
% hold all
% 
% n = 5;
% 
% plotting(A3, A2, P_R_p1_n5+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% legend off
% 
% 
% 
% subplot(2,2,4)
% 
% hold all
% 
% plotting(A3, A2, P_longley_p1+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% figure;
% 
% subplot(2,2,1);
% 
% hold all;
% 
% n = 2;
% 
% plotting(A3, A2, P_R_p2_n2+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% 
% legend off
% title off
% hold all
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% legend off
% 
% 
% subplot(2,2,3)
% 
% hold all
% 
% n = 5;
% 
% plotting(A3, A2, P_R_p2_n5+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% title off
% hold all
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% legend off
% 
% 
% 
% subplot(2,2,4)
% 
% hold all
% 
% plotting(A3, A2, P_longley_p2+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 


% 
% 
% figure;
% 
% subplot(1,2,1);
% 
% hold all;
% 
% n = 2;
% 
% plotting(A3, A2, P_R_p1_n2+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% subplot(1,2,2)
% 
% hold all
% plotting(A3, A2, P_R_p1_n5+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% figure;
% 
% subplot(1,2,1)
% 
% n = 5;
% 
% plotting(A3, A2, P_R_p1_n5+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% subplot(1,2,2)
% 
% hold all
% plotting(A3, A2, P_R_p1_n5+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Power Law n = ', num2str(n)];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% %longley rice
% 
% figure;
% 
% subplot(1,2,1)
% 
% plotting(A3, A2, P_longley_p1+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f1)
% title off
% hold all
% a = ['f = ', num2str(f1/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)
% 
% 
% subplot(1,2,2)
% 
% hold all
% plotting(A3, A2, P_longley_p2+30, XGrid, YGrid, theta_1, phi_1,...
%     theta_pem, phi_pem, theta_roa, phi_roa, theta_b, phi_b, theta_c, phi_c,...
%     x_lower, x_upper, y_lower, y_upper, clim_lower, clim_upper,...
%     xRange, yRange, n, f2)
% a = ['f = ', num2str(f2/(1e6)), ' MHz, Longley-Rice'];
% title(a);
% hold all
% 
% clim auto
% hc = colorbar('EastOutside');
% 
% ylabel(hc,'Received Power [dBm]','FontSize',16)





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

    %figure;

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

function pl = display_pathloss(pm, theta_b, phi_b, tx)

rx = rxsite(Latitude=phi_b,Longitude= theta_b);

pl = pathloss(pm,rx,tx);

end

function R = dist(A2, A3, phi_1,theta_1,r)

Delta_phi = A2-phi_1;

Delta_theta = A3-theta_1;

a = sin(Delta_phi/2)^2+cos(phi_1)*cos(A2)*sin(Delta_theta/2)^2;

c = 2*atan2(sqrt(a),sqrt(1-a));

R = r*c;

end


function [R_b1, R_b2] = break_point(lambda1, h_0, H_0, A2, A3, A4,...
    phi_1,theta_1, avg_terrain_height, r)

R_b1 = ((4*pi/lambda1)*(h_0)*(A4-avg_terrain_height));

if(R_b1 < 0)

    R = dist(A2, A3, phi_1,theta_1,r);

    R_b1 = A4*R/(H_0+A4);

    R_b2 = R_b1;

else

    R = dist(A2, A3, phi_1,theta_1,r);

    R_b1 = h_0*R/(H_0+A4);

    R_b2 = R_b1;

end

end


function SNR = input_SNR(P_T, G_R, G_T, L_p, k, T_ext, T_RX, B_RX,  epsilon_T, eta_T, epsilon_R, eta_R, epsilon_p, F)

num = epsilon_p*P_T*eta_T*epsilon_T*G_T*G_R;

den = k*(T_ext+T_RX/(eta_R*epsilon_R))*B_RX*L_p;

%T_0 = 290;
%
%T_sys = T_0*(10^(0.1*F)-1);
%
%den = k*T_sys*B_RX*L_p;

SNR = 10*log10(abs(num/den));

end


function P_R = received_power(P_T, G_T, G_R, L_p)


P_R = P_T*G_R*G_T./L_p;

end


% function A4 = HAAT(A2, A3, A4, r, reference)
%
% temp_h = zeros(1,length(A2));
%
% for counter = 1:length(A3)
%
%     ref_lat = A2(counter);
%
%     ref_long = A3(counter);
%
%     sum = 0;
%
%     index = 0;
%
%     for counter2 = 1:length(A3)
%
%         Delta_phi = A2(counter2)-ref_lat;
%
%         Delta_theta = A3(counter2)-ref_long;
%
%         a = sin(Delta_phi/2)^2+cos(ref_lat)*cos(A2(counter))*sin(Delta_theta/2)^2;
%
%         c = 2*atan2(sqrt(a),sqrt(1-a));
%
%         dist = r*c;
%
%         if (dist <= reference)
%
%             sum = sum + A4(counter2);
%
%             index = index + 1;
%
%         end
%
%         temp_h(counter) = sum/index;%
%
%     end
%
% end
%
% A4 = A4-temp_h';
%
% end
