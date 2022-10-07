%% Setting output variables

ax = O(1,:);
ay = O(2,:);

Throttle = O(7,:);
Brakes   = O(8,:);

Speed    = O(9,:); 

sliding_energy_lat1 = O(10,:);
sliding_energy_lat2 = O(11,:);
sliding_energy_lat3 = O(12,:);
sliding_energy_lat4 = O(13,:);

sliding_energy_long1 = O(14,:);
sliding_energy_long2 = O(15,:);
sliding_energy_long3 = O(16,:);
sliding_energy_long4 = O(17,:);

sliding_power_lat1 = abs(O(46,:));
sliding_power_lat2 = abs(O(47,:));
sliding_power_lat3 = abs(O(48,:));
sliding_power_lat4 = abs(O(49,:));

sliding_power_long1 = abs(O(50,:));
sliding_power_long2 = abs(O(51,:));
sliding_power_long3 = abs(O(52,:));
sliding_power_long4 = abs(O(53,:));

alpha1 = O(30,:);
alpha2 = O(31,:);
alpha3 = O(32,:);
alpha4 = O(33,:);
kappa1 = O(34,:);
kappa2 = O(35,:);
kappa3 = O(36,:);
kappa4 = O(37,:);

Fz1 = O(18,:);
Fz2 = O(19,:);
Fz3 = O(20,:);
Fz4 = O(21,:);

Fy1 = O(22,:);
Fy2 = O(23,:);
Fy3 = O(24,:);
Fy4 = O(25,:);

YM_fx = O(38,:);
YM_fy = O(39,:);
YM_mz = O(40,:);

stbi = O(41,:);

% eigenvalues = O(42:43,:);


% xModel = Track.xCar - x_star(1,:).*sin(Track.aYaw);
% yModel = Track.yCar + x_star(1,:).*cos(Track.aYaw);

thetaModel = Track.aYaw + x_star(2,:);

xModel = Track.xCar - x_star(1,:).*sin(thetaModel);
yModel = Track.yCar + x_star(1,:).*cos(thetaModel);


%% ploting trajectory

figure; plotbrowser;
plot(Track.xCar, Track.yCar, 'black--', Track.xCarLeft, Track.yCarLeft, 'black', Track.xCarRight, Track.yCarRight, 'black')
% plot(Track.xCarLeft, Track.yCarLeft, 'black', Track.xCarRight, Track.yCarRight, 'black')
title('Optimized Trajectory','Interpreter','latex')
xlabel('X [m]','Interpreter','latex')
ylabel('Y [m]','Interpreter','latex')
grid on
% axis equal
xlim([min(Track.xCar)-50 max(Track.xCar)+50])
daspect([1.5 1.15 10])
hold on
% surf([xModel(:) xModel(:)], [yModel(:) yModel(:)], [x(1:6:end) x(1:6:end)], ...  % Reshape and replicate data
%      'FaceColor', 'none', ...    % Don't bother filling faces with color
%      'EdgeColor', 'interp', ...  % Use interpolated color for edges
%      'LineWidth', 2);            % Make a thicker line
% plot(xModel(1:end-5),yModel(1:end-5))
% plot(xModel(end-5:end),yModel(end-5:end))
plot(xModel,yModel,'LineWidth',2)
hold on
plot(xModel(1),yModel(1),'>','LineWidth',3)
legend('Centerline','Left Track Limit','Right Track Limit','Optimal Trajectory','Start/Finish','Location','best','Interpreter','latex')
grid minor
set(gca,'FontSize',14)
% plot(xModel,yModel)
% view(2);   % Default 2-D view
% colorbar;  % Add a colorbar
% colorbar('westoutside','AxisLocation',rad2deg(x_star(10,:)));

%%
% u
figure; plotbrowser;
plot(Track.sLap,3.6.*Speed,'-','LineWidth',2,'Color',[1 0.2 0])
ylabel('Speed [kmph]','Interpreter','latex')
title('Speed Profile','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
set(gca,'FontSize',14)
grid minor

%%

% ax,ay and r
figure; plotbrowser;
plot(Track.sLap,ay./9.806,'LineWidth',2)
ylabel('Lateral Acceleration [g]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Lateral Acceleration','Interpreter','latex')
set(gca,'FontSize',14)
grid minor

figure; plotbrowser;
plot(Track.sLap,ax./9.806,'LineWidth',2)
ylabel('Longitudinal Acceleration [g]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Longitudinal Acceleration','Interpreter','latex')
set(gca,'FontSize',14)
grid minor

figure; plotbrowser;
plot(Track.sLap,sqrt((ay./9.806).^2+(ax./9.806).^2),'LineWidth',2,'Color',[1 0 0])
ylabel('Combined Acceleration [g]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Combined Acceleration','Interpreter','latex')
set(gca,'FontSize',14)
grid minor
%%

% figure; plotbrowser;
plot(ay./9.806,ax./9.806,'.','MarkerSize',12)
title('g-g Diagram','Interpreter','latex')
ylabel('Long G [g]','Interpreter','latex')
xlabel('Lat G [g]','Interpreter','latex')
ylim([-2.5 2.5])
xlim([-2.5 2.5])
set(gca,'FontSize',14)
grid minor

%% yaw rate

% figure; plotbrowser;
% plot(Track.sLap,rad2deg(x_star(5,:)))
% title('Yaw Velocity profile over a lap')
% set(gca,'FontSize',14)
% grid minor

%% wheel speeds
% figure; plotbrowser;
% subplot(1,2,1)
% plot(Track.sLap,3.6.*(x_star(6,:).*Vehicle.tire_1.radius) )
% hold on
% plot(Track.sLap,3.6.*(x_star(7,:).*Vehicle.tire_2.radius) )
% ylabel('kmph')
% title('front - wheel speed')
% legend('FL','FR')
% subplot(1,2,2)
% plot(Track.sLap,3.6.*(x_star(8,:).*Vehicle.tire_3.radius) )
% hold on
% plot(Track.sLap,3.6.*(x_star(9,:).*Vehicle.tire_4.radius) )
% ylabel('kmph')
% title('rear - wheel speed')
% legend('RL','RR')

%% Steering

figure; plotbrowser;
plot(Track.sLap,rad2deg(x_star(10,:)),'LineWidth',2)
title('Optimal Steering Angle Profile','Interpreter','latex')
ylabel('Steering Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
set(gca,'FontSize',14)
grid minor

%% Tau

figure; plotbrowser;
subplot(2,1,1)
plot(Track.sLap,100.*Throttle,'LineWidth',2)
ylabel('Percentage [\%]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Throttle Position Profile Comparison','Interpreter','latex')
set(gca,'FontSize',14)
grid minor
subplot(2,1,2)
plot(Track.sLap,-5000.*Brakes,'LineWidth',2)
ylabel('Percentage [\%]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Brake Torque Profile Comparison','Interpreter','latex')
set(gca,'FontSize',14)
grid minor

%% control rates

figure; plotbrowser;
subplot(2,1,1)
plot(Track.sLap,(u_star(1,:)))
title('delta rate')
subplot(2,1,2)
plot(Track.sLap,(u_star(2,:)))
title('tau rate')

%% Sliding energy lateral

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,sliding_energy_lat1/1e3,'--','Color',[1 0 0],'LineWidth',2)
title('Sliding Energy Lateral FL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,2)
plot(Track.sLap,sliding_energy_lat2/1e3,'--','Color',[0 1 0],'LineWidth',2)
title('Sliding Energy Lateral FR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,3)
plot(Track.sLap,sliding_energy_lat3/1e3,'--','Color',[0.3 0.3 1],'LineWidth',2)
title('Sliding Energy Lateral RL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,4)
plot(Track.sLap,sliding_energy_lat4/1e3,'--','Color',[1 1 0.0],'LineWidth',2)
title('Sliding Energy Lateral RR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%% Sliding energy longitudinal

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,sliding_energy_long1/1e3,'--','Color',[1 0 0],'LineWidth',2)
title('Sliding Energy Longitudinal FL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,2)
plot(Track.sLap,sliding_energy_long2/1e3,'--','Color',[0 1 0],'LineWidth',2)
title('Sliding Energy Longitudinal FR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,3)
plot(Track.sLap,sliding_energy_long3/1e3,'--','Color',[0.3 0.3 1],'LineWidth',2)
title('Sliding Energy Longitudinal RL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,4)
plot(Track.sLap,sliding_energy_long4/1e3,'--','Color',[1 1 0.0],'LineWidth',2)
title('Sliding Energy Longitudinal RR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%% Sliding power lateral

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,sliding_power_lat1/1e3,'--','Color',[1 0 0],'LineWidth',2)
title('Sliding power Lateral FL','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,2)
plot(Track.sLap,sliding_power_lat2/1e3,'--','Color',[0 1 0],'LineWidth',2)
title('Sliding power Lateral FR','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,3)
plot(Track.sLap,sliding_power_lat3/1e3,'--','Color',[0.3 0.3 1],'LineWidth',2)
title('Sliding power Lateral RL','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,4)
plot(Track.sLap,sliding_power_lat4/1e3,'--','Color',[1 1 0.0],'LineWidth',2)
title('Sliding power Lateral RR','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%% Sliding power longitudinal

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,sliding_power_long1/1e3,'--','Color',[1 0 0],'LineWidth',2)
title('Sliding power Longitudinal FL','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,2)
plot(Track.sLap,sliding_power_long2/1e3,'--','Color',[0 1 0],'LineWidth',2)
title('Sliding power Longitudinal FR','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,3)
plot(Track.sLap,sliding_power_long3/1e3,'--','Color',[0.3 0.3 1],'LineWidth',2)
title('Sliding power Longitudinal RL','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%%
subplot(2,2,4)
plot(Track.sLap,sliding_power_long4/1e3,'--','Color',[1 1 0.0],'LineWidth',2)
title('Sliding power Longitudinal RR','Interpreter','latex')
ylabel('power [kW]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%% lateral forces

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,Fy1,'Color',[1 0 0],'LineWidth',2)
title('Lateral Force  FL')
ylabel('Lateral Force [N]')
xlabel('Distance [m]')
grid minor
set(gca,'FontSize',13)
subplot(2,2,2)
plot(Track.sLap,Fy2,'Color',[0 1 0],'LineWidth',2)
title('Lateral Force  FR')
ylabel('Lateral Force [N]')
xlabel('Distance [m]')
grid minor
set(gca,'FontSize',13)
subplot(2,2,3)
plot(Track.sLap,Fy3,'Color',[0.3 0.3 1],'LineWidth',2)
title('Lateral Force  RL')
ylabel('Lateral Force [N]')
xlabel('Distance [m]')
grid minor
set(gca,'FontSize',13)
subplot(2,2,4)
plot(Track.sLap,Fy4,'Color',[1 1 0.0],'LineWidth',2)
title('Lateral Force  RR')
ylabel('Lateral Force [N]')
xlabel('Distance [m]')
grid minor
set(gca,'FontSize',13)

%% Fz

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,Fz1,'Color',[1 0 0],'LineWidth',2)
title('Tire Normal Load - FL','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
% ylim([0*min(Fz1) max(Fz2)])
ylim([0*min(Fz1) 8e3])
grid minor
set(gca,'FontSize',13)
subplot(2,2,2)
plot(Track.sLap,Fz2,'Color',[0.2 0.8 0],'LineWidth',2)
title('Tire Normal Load - FR','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
ylim([0*min(Fz1) 8e3])
grid minor
set(gca,'FontSize',13)
subplot(2,2,3)
plot(Track.sLap,Fz3,'Color',[0.3 0.3 1],'LineWidth',2)
title('Tire Normal Load - RL','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
ylim([0*min(Fz1) 8e3])
grid minor
set(gca,'FontSize',13)
subplot(2,2,4)
plot(Track.sLap,Fz4,'Color',[1 0.8 0.0],'LineWidth',2)
title('Tire Normal Load - RR','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
ylim([0*min(Fz1) 8e3])
grid minor
set(gca,'FontSize',13)

%% stability index

figure;plotbrowser;
plot(Track.sLap,stbi)
title('Static Margin')
ylabel('Static Margin [-]')
xlabel('Distance [m]')
yline(0)
grid minor
set(gca,'FontSize',15)

%% ym
figure; plotbrowser;
plot(Track.sLap,YM_fy,'LineWidth',2)
hold on
plot(Track.sLap,YM_fx,'LineWidth',2)
plot(Track.sLap,YM_mz,'LineWidth',2)
plot(Track.sLap,YM_fy+YM_fx+YM_mz,'--','LineWidth',2)
ylabel('Yaw Moment [Nm]')
xlabel('Distance [m]')
legend('Yaw Moment From Fy','Yaw Moment From Fx','Yaw Moment From Mz','Total')
title('3 Different Yaw Moment profile over a Lap')
set(gca,'FontSize',14)
grid minor

%%
figure; plotbrowser;
plot(Track.sLap,YM_fy+YM_fx+YM_mz,'-','LineWidth',2)
ylabel('Yaw Moment [Nm]')
xlabel('Distance [m]')
title('Yaw Moment profile over a Lap')
set(gca,'FontSize',14)

%% ym vs lat g
figure; plotbrowser;
plot(ay/9.806,YM_fy+YM_fx+YM_mz,'o','LineWidth',2)
ylabel('Yaw Moment [Nm]')
xlabel('Lateral Acceleration [g]')
title('Yaw Moment profile over a Lap')
set(gca,'FontSize',14)

%% alphas

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,57.29*alpha1,'--','Color',[1 0 0],'LineWidth',2)
title('Slip Angle FL','Interpreter','latex')
ylabel('Slip Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

subplot(2,2,2)
plot(Track.sLap,57.29*alpha2,'--','Color',[0 1 0],'LineWidth',2)
title('Slip Angle FR','Interpreter','latex')
ylabel('Slip Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

subplot(2,2,3)
plot(Track.sLap,57.29*alpha3,'--','Color',[0.3 0.3 1],'LineWidth',2)
title('Slip Angle RL','Interpreter','latex')
ylabel('Slip Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

subplot(2,2,4)
plot(Track.sLap,57.29*alpha4,'--','Color',[1 1 0.0],'LineWidth',2)
title('Slip Angle RR','Interpreter','latex')
ylabel('Slip Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)