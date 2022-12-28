%testing the vehicle model
clc
clear

x0 = zeros(30,1);

vx0 = 160/3.6;
x0(3,1) = -0.03;
x0(4,1) = 0.0;
x0(5,1) = 0;
x0(7,1) = -0.01;
x0(8,1) = -0.01;
x0(9,1) = -0.01;
x0(10,1) = -0.01;
x0(15,1) = vx0; % longitudinal speed
x0(25,1) = vx0/0.32; % longitudinal speed
x0(26,1) = vx0/0.32; % longitudinal speed
x0(27,1) = vx0/0.32; % longitudinal speed
x0(28,1) = vx0/0.32; % longitudinal speed
x0(29,1) = deg2rad(20); % steer
x0(30,1) = 0; % tau

tspan = 0:0.01:3;

u = [0 0]';

p = 0.5*ones(4,1);
p(3) = 7e3;

vehicle = Model.load_vehicle_parameters();


% [tout, yout] = ode45(@(t,x) DynamicsOpenLoop(t,x,u,vehicle),tspan,x0);
[tout, yout] = ode45(@(t,x) Controller.fnDynamics14DOFVehicleOLSimple(t,x,u,p,vehicle),tspan,x0);

% u_star = zeros(4,length(tout));

% [~,~,~,O] = fnDynamicsVehicle_linear_tires([],yout',u_star,vehicle,p);

% dx = zeros(3,length(tout));

% [dx,~,~] = Controller.fnDynamicsVehicle([],yout',u,vehicle);

%% ploting

yaw_rate = yout(:,20);
long_speed = 3.6.*yout(:,15);

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
