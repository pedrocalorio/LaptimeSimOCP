%testing the vehicle model
clc
clear

x0 = zeros(9,1);
vx0 = 150/3.6;
x0(1,1) = vx0; % longitudinal speed
x0(4,1) = vx0/0.33; % longitudinal speed
x0(5,1) = vx0/0.33; % longitudinal speed
x0(6,1) = vx0/0.33; % longitudinal speed
x0(7,1) = vx0/0.33; % longitudinal speed
x0(8,1) = deg2rad(2); % steer
x0(9,1) = 0; % tau

tspan = 0:0.01:2;

u = [0 0 0 0]';

p = 0.5*ones(4,1);

vehicle = Model.load_vehicle_parameters();

% [tout, yout] = ode45(@(t,x) DynamicsOpenLoop(t,x,u,vehicle),tspan,x0);
[tout, yout] = ode15s(@(t,x) fnDynamicsVehicle_linear_tires(t,x,u,vehicle,p),tspan,x0);

u_star = zeros(4,length(tout));

[~,~,~,O] = fnDynamicsVehicle_linear_tires([],yout',u_star,vehicle,p);

% dx = zeros(3,length(tout));

% [dx,~,~] = Controller.fnDynamicsVehicle([],yout',u,vehicle);

%% ploting

yaw_rate = yout(:,3);
long_speed = 3.6.*yout(:,1);

% figure;
% plot(tout,rad2deg(yaw_rate))
% title('Yaw Rate')

figure;
plot(tout,long_speed)
title('Longitudinal Speed')
hold on

figure;
plot(tout,yaw_rate)
hold on

% %% pack unpack
% 
% [z,pack] = UpackDecVar(yout',u_star);
% 
% [x_nw,u_nw] = unPackDecVar(z,pack);
