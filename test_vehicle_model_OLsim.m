%testing the vehicle model
clc
clear
close all

x0 = zeros(38,1);
vx0 = 150/3.6;

x0(3,1) = 0.01;
x0(4,1) = 0.0;
x0(5,1) = 0.02/57.3;
x0(7,1) = 0.01;
x0(8,1) = 0.01;
x0(9,1) = 0.01;
x0(10,1) = 0.01;

x0(15,1) = vx0; % longitudinal speed

x0(25,1) = vx0/0.32; % longitudinal speed
x0(26,1) = vx0/0.32; % longitudinal speed
x0(27,1) = vx0/0.32; % longitudinal speed
x0(28,1) = vx0/0.32; % longitudinal speed

x0(29,1) = deg2rad(30); % steer
x0(30,1) = 0*0.06; % tau

tspan = 0:0.01:3;

u = [0 0]';

p = 0.5*ones(4,1);
p(3) = 7e3;

vehicle = Model.load_vehicle_parameters();

m = vehicle.chassis.Ms; % Sprung mass (kg)
g = 9.81;
wheelbase   = vehicle.chassis.wheelbase;
a           = wheelbase * (1-p(1));
b           = wheelbase * p(1);
muf=55;    %front unsprung mass (kg)
mur=55;    %rear unsprung mass (kg)

kslf=   160e3; %front left suspension stiffness (N/m)
ksrf=   160e3;
kslr=   95e3;  %rear left suspension stiffness (N/m)
ksrr=   95e3;
ktf=300e3;  %front tire stiffness (N/m)
ktr=300e3;  %rear tire stiffness (N/m)

xtirf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
xtilf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
xtilr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
xtirr=((m*g*a)/(2*(a+b))+mur*g)/ktr;

xsirf=(m*g*b)/(2*(a+b)*ksrf);
xsilf=(m*g*b)/(2*(a+b)*kslf);
xsilr=(m*g*a)/(2*(a+b)*kslr);
xsirr=(m*g*a)/(2*(a+b)*ksrr);

x0(31,1) = 15.37/1000;
x0(32,1) = 15.37/1000;
x0(33,1) = 27.26/1000;
x0(34,1) = 27.26/1000;
x0(35,1) = 7.71/1000;
x0(36,1) = 7.71/1000;
x0(37,1) = 8.176/1000;
x0(38,1) = 8.176/1000;

% [tout, yout] = ode45(@(t,x) DynamicsOpenLoop(t,x,u,vehicle),tspan,x0);
[tout, yout] = ode15s(@(t,x) Controller.fnDynamics14DOFVehicleOL(t,x,u,p,vehicle),tspan,x0);

% u_star = zeros(4,length(tout));

% [~,~,~,O] = fnDynamicsVehicle_linear_tires([],yout',u_star,vehicle,p);

% dx = zeros(3,length(tout));

% [dx,~,~] = Controller.fnDynamicsVehicle([],yout',u,vehicle);

%% ploting

yaw_rate    = yout(:,20);
long_speed  = 3.6.*yout(:,15);

roll_angle  = rad2deg(yout(:,4));
pitch_angle = rad2deg(yout(:,5));
psi_angle   = rad2deg(yout(:,6));

xsfl = 1000*yout(:,31);
xsfr = 1000*yout(:,32);
xsrl = 1000*yout(:,33);
xsrr = 1000*yout(:,34);

xtfl = 1000*yout(:,35);
xtfr = 1000*yout(:,36);
xtrl = 1000*yout(:,37);
xtrr = 1000*yout(:,38);

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

figure;
plot(tout,roll_angle)
title('Roll Angle')
hold on

figure;
plot(tout,pitch_angle)
title('Pitch Angle')
hold on

figure;
plot(tout,psi_angle)
title('Psi Angle')
hold on


figure;
plot(yout(:,1),yout(:,2))
title('Vehicle Displacement')
hold on

figure;
subplot(2,2,1)
plot(tout,xsfl)
subplot(2,2,2)
plot(tout,xsfr)
subplot(2,2,3)
plot(tout,xsrl)
subplot(2,2,4)
plot(tout,xsrr)

figure;
subplot(2,2,1)
plot(tout,xtfl)
subplot(2,2,2)
plot(tout,xtfr)
subplot(2,2,3)
plot(tout,xtrl)
subplot(2,2,4)
plot(tout,xtrr)



% %% pack unpack
% 
% [z,pack] = UpackDecVar(yout',u_star);
% 
% [x_nw,u_nw] = unPackDecVar(z,pack);
