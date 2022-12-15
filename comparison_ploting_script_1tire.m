baseline = load("SimResults/tire energy saving study - fl tires (combined)/budapest_baseline/sim.mat");

saving_energy_10 = load("SimResults/tire energy saving study - fl tires (combined)/budapest_90p_fl_tire_comb_energy/sim.mat");

saving_energy_20 = load("SimResults/tire energy saving study - fl tires (combined)/budapest_80p_fl_tire_comb_energy/sim.mat");

saving_energy_30 = load("SimResults/tire energy saving study - fl tires (combined)/budapest_70p_fl_tire_comb_energy/sim.mat");

saving_energy_40 = load("SimResults/tire energy saving study - fl tires (combined)/budapest_60p_fl_tire_comb_energy/sim.mat");


Track = baseline.Track;

%%
figure; plotbrowser;
energy_array = [baseline.KPIs.combined_sliding_energy_1 saving_energy_10.KPIs.combined_sliding_energy_1 ...
    saving_energy_20.KPIs.combined_sliding_energy_1 saving_energy_30.KPIs.combined_sliding_energy_1...
    saving_energy_40.KPIs.combined_sliding_energy_1];
laptime_array = [baseline.KPIs.laptime saving_energy_10.KPIs.laptime ...
    saving_energy_20.KPIs.laptime saving_energy_30.KPIs.laptime...
    saving_energy_40.KPIs.laptime];
plot(energy_array/1e3,laptime_array,'-o')
xlabel('Combined Sliding Energy (FL tire) [kJ]')
ylabel('Laptime [s]')
title('Laptime vs Combined Sliding Energy Front Left Tire')

%%

% speed
figure; plotbrowser;
plot(Track.sLap,3.6.*baseline.O(9,:),'-','LineWidth',1,'Color',[1 0.2 0])
hold on
plot(Track.sLap,3.6.*saving_energy_10.O(9,:),'-','LineWidth',1,'Color','b')
plot(Track.sLap,3.6.*saving_energy_20.O(9,:),'-','LineWidth',1,'Color','g')
plot(Track.sLap,3.6.*saving_energy_30.O(9,:),'-','LineWidth',1)
plot(Track.sLap,3.6.*saving_energy_40.O(9,:),'-','LineWidth',1)
ylabel('Speed [kmph]','Interpreter','latex')
title('Speed Profile','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
legend('Baseline','10% less tire energy','20% less tire energy',...
    '30% less tire energy','40% less tire energy',...
    'location','best')
set(gca,'FontSize',14)
grid minor
%%
% ax,ay and r
figure; plotbrowser;
plot(Track.sLap,baseline.O(2,:)./9.806,'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(2,:)./9.806,'LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(2,:)./9.806,'LineWidth',2,'Color','g')
ylabel('Lateral Acceleration [g]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Lateral Acceleration','Interpreter','latex')
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
set(gca,'FontSize',14)
grid minor

figure; plotbrowser;
plot(Track.sLap,baseline.O(1,:)./9.806,'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(1,:)./9.806,'LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(1,:)./9.806,'LineWidth',2,'Color','g')
ylabel('Longitudinal Acceleration [g]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Longitudinal Acceleration','Interpreter','latex')
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
set(gca,'FontSize',14)
grid minor

figure; plotbrowser;
plot(baseline.O(2,:)./9.806,baseline.O(1,:)./9.806,'.','MarkerSize',12,'Color',[1 0.2 0])
hold on
plot(saving_energy_15.O(2,:)./9.806,saving_energy_15.O(1,:)./9.806,'.','MarkerSize',12,'Color','b')
plot(saving_energy_30.O(2,:)./9.806,saving_energy_30.O(1,:)./9.806,'.','MarkerSize',12,'Color','g')
title('g-g Diagram','Interpreter','latex')
ylabel('Long G [g]','Interpreter','latex')
xlabel('Lat G [g]','Interpreter','latex')
ylim([-2.5 2.5])
xlim([-2.5 2.5])
set(gca,'FontSize',14)
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
grid minor

%% yaw rate

figure; plotbrowser;
plot(Track.sLap,rad2deg(baseline.x_star(5,:)),'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,rad2deg(saving_energy_15.x_star(5,:)),'LineWidth',2,'Color','b')
plot(Track.sLap,rad2deg(saving_energy_30.x_star(5,:)),'LineWidth',2,'Color','g')
ylabel('Yaw Rate [deg/s]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Yaw Velocity profile over a lap')
set(gca,'FontSize',14)
grid minor

%% Steering

figure; plotbrowser;
plot(Track.sLap,rad2deg(baseline.x_star(10,:)),'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,rad2deg(saving_energy_15.x_star(10,:)),'LineWidth',2,'Color','b')
plot(Track.sLap,rad2deg(saving_energy_30.x_star(10,:)),'LineWidth',2,'Color','g')
title('Optimal Steering Angle Profile','Interpreter','latex')
ylabel('Steering Angle [deg]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
set(gca,'FontSize',14)
grid minor

%% Tau

figure; plotbrowser;
subplot(2,1,1)
plot(Track.sLap,100.*baseline.O(7,:),'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,100.*saving_energy_15.O(7,:),'LineWidth',2,'Color','b')
plot(Track.sLap,100.*saving_energy_30.O(7,:),'LineWidth',2,'Color','g')
ylabel('Percentage [\%]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Throttle Position Profile Comparison','Interpreter','latex')
set(gca,'FontSize',14)
grid minor
subplot(2,1,2)
plot(Track.sLap,-5000.*baseline.O(8,:),'LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,-5000.*saving_energy_15.O(8,:),'LineWidth',2,'Color','b')
plot(Track.sLap,-5000.*saving_energy_30.O(8,:),'LineWidth',2,'Color','g')
ylabel('Torque [Nm]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
title('Brake Torque Profile Comparison','Interpreter','latex')
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
set(gca,'FontSize',14)
grid minor

%% Sliding energy lateral

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,baseline.O(10,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(10,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(10,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Lateral FL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,2)
plot(Track.sLap,baseline.O(11,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(11,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(11,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Lateral FR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,3)
plot(Track.sLap,baseline.O(12,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(12,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(12,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Lateral RL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,4)
plot(Track.sLap,baseline.O(13,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(13,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(13,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Lateral RR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%% Sliding energy long

figure;plotbrowser;
subplot(2,2,1)
plot(Track.sLap,baseline.O(14,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(14,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(14,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Longitudinal FL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,2)
plot(Track.sLap,baseline.O(15,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(15,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(15,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Longitudinal FR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,3)
plot(Track.sLap,baseline.O(16,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(16,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(16,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Longitudinal RL','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)
%
subplot(2,2,4)
plot(Track.sLap,baseline.O(17,:)/1e3,'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(17,:)/1e3,'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(17,:)/1e3,'-','LineWidth',2,'Color','g')
title('Sliding Energy Longitudinal RR','Interpreter','latex')
ylabel('Energy [kJ]','Interpreter','latex')
xlabel('Distance [m]','Interpreter','latex')
grid minor
set(gca,'FontSize',13)

%%
figure; plotbrowser;
plot(Track.sLap,baseline.O(38,:)+baseline.O(39,:)+baseline.O(40,:),'-','LineWidth',2,'Color',[1 0.2 0])
hold on
plot(Track.sLap,saving_energy_15.O(38,:)+saving_energy_15.O(39,:)+saving_energy_15.O(40,:),'-','LineWidth',2,'Color','b')
plot(Track.sLap,saving_energy_30.O(38,:)+saving_energy_30.O(39,:)+saving_energy_30.O(40,:),'-','LineWidth',2,'Color','g')
ylabel('Yaw Moment [Nm]')
xlabel('Distance [m]')
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
title('Yaw Moment profile over a Lap')
grid minor
set(gca,'FontSize',14)

%% ym vs lat g
figure; plotbrowser;
plot(baseline.O(2,:)/9.806,baseline.O(38,:)+baseline.O(39,:)+baseline.O(40,:),'.','MarkerSize',16,'Color',[1 0.2 0])
hold on
plot(saving_energy_15.O(2,:)/9.806,saving_energy_15.O(38,:)+saving_energy_15.O(39,:)+saving_energy_15.O(40,:),'.','MarkerSize',16,'Color','b')
plot(saving_energy_30.O(2,:)/9.806,saving_energy_30.O(38,:)+saving_energy_30.O(39,:)+saving_energy_30.O(40,:),'.','MarkerSize',16,'Color','g')
ylabel('Yaw Moment [Nm]')
xlabel('Lateral Acceleration [g]')
xlim([-3.5 3.5])
legend('Baseline','15% less tire energy','30% less tire energy','location','best')
title('Yaw Moment profile over a Lap')
set(gca,'FontSize',14)
grid minor