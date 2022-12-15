function KPI = metrics_output(O,problem,soln)

Track = problem.dsSystem.td;
Vehicle = problem.dsSystem.vd;





n = soln.state(1,:);
zeta = soln.state(2,:);
vx = soln.state(3,:);
vy = soln.state(4,:);
% r = soln.grid.state(5,:);

curv = interp1(Track.distance,Track.curv,Track.sLap);
Sf = (1 - n.*curv)./(vx.*cos(zeta)-vy.*sin(zeta));
KPI.laptime = trapz(Track.sLap,Sf);

% steering integral
KPI.steering_Integral = (trapz(Track.sLap,abs(rad2deg(soln.state(10,:))) ) );
KPI.throttle_Integral = trapz(Track.sLap,O(7,:));
KPI.brakes_Integral   = trapz(Track.sLap,O(8,:));

%%%%%%%%%%%%%%%
% sliding energy
%%%%%%%%%%%%%%%

% lateral
KPI.lateral_sliding_energy_1 = trapz(O(1,:),O(10,:))/1e3;
KPI.lateral_sliding_energy_2 = trapz(O(1,:),O(11,:))/1e3;
KPI.lateral_sliding_energy_3 = trapz(O(1,:),O(12,:))/1e3;
KPI.lateral_sliding_energy_4 = trapz(O(1,:),O(13,:))/1e3;

KPI.lateral_sliding_energy_1_dist = O(57,1);
KPI.lateral_sliding_energy_2_dist = O(58,1);
KPI.lateral_sliding_energy_3_dist = O(59,1);
KPI.lateral_sliding_energy_4_dist = O(60,1);

KPI.lateral_sliding_energy_front = KPI.lateral_sliding_energy_1+KPI.lateral_sliding_energy_2;
KPI.lateral_sliding_energy_rear = KPI.lateral_sliding_energy_3+KPI.lateral_sliding_energy_4;
KPI.lateral_sliding_energy_distribution = KPI.lateral_sliding_energy_front/(KPI.lateral_sliding_energy_front+KPI.lateral_sliding_energy_rear);

% longitudinal
KPI.long_sliding_energy_1 = trapz(O(1,:),O(14,:))/1e3;
KPI.long_sliding_energy_2 = trapz(O(1,:),O(15,:))/1e3;
KPI.long_sliding_energy_3 = trapz(O(1,:),O(16,:))/1e3;
KPI.long_sliding_energy_4 = trapz(O(1,:),O(17,:))/1e3;

KPI.long_sliding_energy_1_dist = O(61,1);
KPI.long_sliding_energy_2_dist = O(62,1);
KPI.long_sliding_energy_3_dist = O(63,1);
KPI.long_sliding_energy_4_dist = O(64,1);


KPI.long_sliding_energy_front = KPI.long_sliding_energy_1+KPI.long_sliding_energy_2;
KPI.long_sliding_energy_rear = KPI.long_sliding_energy_3+KPI.long_sliding_energy_4;
KPI.long_sliding_energy_distribution = KPI.long_sliding_energy_front/(KPI.long_sliding_energy_front+KPI.long_sliding_energy_rear);

% total
KPI.combined_sliding_energy_1 = KPI.lateral_sliding_energy_1+KPI.long_sliding_energy_1;
KPI.combined_sliding_energy_2 = KPI.lateral_sliding_energy_2+KPI.long_sliding_energy_2;
KPI.combined_sliding_energy_3 = KPI.lateral_sliding_energy_3+KPI.long_sliding_energy_3;
KPI.combined_sliding_energy_4 = KPI.lateral_sliding_energy_4+KPI.long_sliding_energy_4;

KPI.combined_sliding_energy_1_dist = KPI.lateral_sliding_energy_1_dist+KPI.long_sliding_energy_1_dist;
KPI.combined_sliding_energy_2_dist = KPI.lateral_sliding_energy_2_dist+KPI.long_sliding_energy_2_dist;
KPI.combined_sliding_energy_3_dist = KPI.lateral_sliding_energy_3_dist+KPI.long_sliding_energy_3_dist;
KPI.combined_sliding_energy_4_dist = KPI.lateral_sliding_energy_4_dist+KPI.long_sliding_energy_4_dist;

KPI.combined_sliding_energy_front = KPI.combined_sliding_energy_1+KPI.combined_sliding_energy_2;
KPI.combined_sliding_energy_rear = KPI.combined_sliding_energy_3+KPI.combined_sliding_energy_4;
KPI.combined_sliding_energy_total = KPI.combined_sliding_energy_front+KPI.combined_sliding_energy_rear;
% KPI.combined_sliding_energy_distribution = KPI.combined_sliding_energy_front/(KPI.combined_sliding_energy_front+KPI.combined_sliding_energy_rear);


KPI.combined_sliding_energy_front_dist = KPI.combined_sliding_energy_1_dist+KPI.combined_sliding_energy_2_dist;
KPI.combined_sliding_energy_rear_dist  = KPI.combined_sliding_energy_3_dist+KPI.combined_sliding_energy_4_dist;
KPI.combined_sliding_energy_total_dist = KPI.combined_sliding_energy_front_dist+KPI.combined_sliding_energy_rear_dist;

KPI.lateral_tire_efficiency_1 = trapz(O(1,:),O(69,:));
KPI.lateral_tire_efficiency_2 = trapz(O(1,:),O(70,:));
KPI.lateral_tire_efficiency_3 = trapz(O(1,:),O(71,:));
KPI.lateral_tire_efficiency_4 = trapz(O(1,:),O(72,:));

KPI.long_tire_efficiency_1 = trapz(O(1,:),O(65,:));
KPI.long_tire_efficiency_2 = trapz(O(1,:),O(66,:));
KPI.long_tire_efficiency_3 = trapz(O(1,:),O(67,:));
KPI.long_tire_efficiency_4 = trapz(O(1,:),O(68,:));


end