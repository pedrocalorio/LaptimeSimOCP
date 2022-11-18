function KPI = metrics_output(O,problem,soln)

Track = problem.dsSystem.td;
Vehicle = problem.dsSystem.vd;

KPI.laptime = soln.obj_func;

n = soln.state(1,:);
zeta = soln.state(2,:);
vx = soln.state(3,:);
vy = soln.state(4,:);
% r = soln.grid.state(5,:);


% steering integral
KPI.steering_Integral = (trapz(abs(rad2deg(soln.state(10,:))) ) );
KPI.throttle_Integral = trapz(O(7,:));
KPI.brakes_Integral   = trapz(O(8,:));

%%%%%%%%%%%%%%%
% sliding energy
%%%%%%%%%%%%%%%

% lateral
KPI.lateral_sliding_energy_1 = trapz(O(45,:),O(10,:))/1e3;
KPI.lateral_sliding_energy_2 = trapz(O(45,:),O(11,:))/1e3;
KPI.lateral_sliding_energy_3 = trapz(O(45,:),O(12,:))/1e3;
KPI.lateral_sliding_energy_4 = trapz(O(45,:),O(13,:))/1e3;

KPI.lateral_sliding_energy_front = KPI.lateral_sliding_energy_1+KPI.lateral_sliding_energy_2;
KPI.lateral_sliding_energy_rear = KPI.lateral_sliding_energy_3+KPI.lateral_sliding_energy_4;
KPI.lateral_sliding_energy_distribution = KPI.lateral_sliding_energy_front/(KPI.lateral_sliding_energy_front+KPI.lateral_sliding_energy_rear);

% longitudinal
KPI.long_sliding_energy_1 = trapz(O(45,:),O(14,:))/1e3;
KPI.long_sliding_energy_2 = trapz(O(45,:),O(15,:))/1e3;
KPI.long_sliding_energy_3 = trapz(O(45,:),O(16,:))/1e3;
KPI.long_sliding_energy_4 = trapz(O(45,:),O(17,:))/1e3;

KPI.long_sliding_energy_front = KPI.long_sliding_energy_1+KPI.long_sliding_energy_2;
KPI.long_sliding_energy_rear = KPI.long_sliding_energy_3+KPI.long_sliding_energy_4;
KPI.long_sliding_energy_distribution = KPI.long_sliding_energy_front/(KPI.long_sliding_energy_front+KPI.long_sliding_energy_rear);

% total
KPI.total_sliding_energy_1 = KPI.lateral_sliding_energy_1+KPI.long_sliding_energy_1;
KPI.total_sliding_energy_2 = KPI.lateral_sliding_energy_2+KPI.long_sliding_energy_2;
KPI.total_sliding_energy_3 = KPI.lateral_sliding_energy_3+KPI.long_sliding_energy_3;
KPI.total_sliding_energy_4 = KPI.lateral_sliding_energy_4+KPI.long_sliding_energy_4;

KPI.total_sliding_energy_front = KPI.total_sliding_energy_1+KPI.total_sliding_energy_2;
KPI.total_sliding_energy_rear = KPI.total_sliding_energy_3+KPI.total_sliding_energy_4;
KPI.total_sliding_energy_distribution = KPI.total_sliding_energy_front/(KPI.total_sliding_energy_front+KPI.total_sliding_energy_rear);



KPI.Sf = (1 - n.*Track.curv)./(vx.*cos(zeta)-vy.*sin(zeta));


end