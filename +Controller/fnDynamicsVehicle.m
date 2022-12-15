function [dx, g_ineq, q_eq,saving_constraints] = fnDynamicsVehicle(x,u,p,Vehicle,Track)

    addpath('C:/dev/libraries/casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*

    kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

    % Distance Step
    ds = (Track.sLap(end) - Track.sLap(1)) / (length(x(1,:)) - 1);
    
    % Gravity constant
    g = 9.81;
    
    %% System states
    
    n       = x(1,:);
    xi      = x(2,:);
    vx      = x(3,:); % velocity
    vy      = x(4,:); % body-slip angle
    r       = x(5,:); % yaw rate
    omega1  = x(6,:);
    omega2  = x(7,:);
    omega3  = x(8,:);
    omega4  = x(9,:);
    steer   = x(10,:); 
    tau     = x(11,:); 

    
    %% Control inputs

    steerRate   = u(1,:);
    tauRate     = u(2,:);
    LongLT_norm = u(3,:);
    LatLT_norm  = u(4,:);    

    %% Calculating the cumulative laptime

    Sf = (1 - n.*kappa)./(vx.*cos(xi)-vy.*sin(xi));
    cumLaptime = ds*Sf*ones(length(n),1);
    
    %% Vehicle Parameters 
    
    % Aero
    liftCoeff   = Vehicle.aero.liftCoeff;
    dragCoeff   = Vehicle.aero.dragCoeff;
    frontalArea = Vehicle.aero.frontalArea;
    airDens     = Vehicle.aero.airDensity;
    aeroBalance = p(2);
    
    % Susp
    LatLTD  = p(3);
    
    % Steering
    steringRatio = Vehicle.steeringratio;
    
    % Chassis
    m           = Vehicle.chassis.Ms;
    wheelbase   = Vehicle.chassis.wheelbase;
    weight_dist = p(1);
    a           = wheelbase * (1-weight_dist);
    b           = wheelbase * weight_dist;
    frontTrack  = Vehicle.chassis.frontTrack;
    rearTrack   = Vehicle.chassis.rearTrack;
    Iz          = Vehicle.chassis.Iz;
    CoG         = Vehicle.chassis.CoG;
    
    % Wheel Inertias
    Iwheel_1 = Vehicle.tire_1.inertia;
    Iwheel_2 = Vehicle.tire_2.inertia;  
    Iwheel_3 = Vehicle.tire_3.inertia;
    Iwheel_4 = Vehicle.tire_4.inertia;  
    
    Rl1 = Vehicle.tire_1.radius;
    Rl2 = Vehicle.tire_2.radius;
    Rl3 = Vehicle.tire_3.radius;
    Rl4 = Vehicle.tire_4.radius;
    
    MF_f = Vehicle.tire_1.MF;
    MF_r = Vehicle.tire_3.MF;        
    
    
    %% Calculations
    
    % Steering angles @ tires
    delta = (steer/steringRatio) ; %radians  
        
    % (un)normalized longitudinal weight transfer
    LongLT = (LongLT_norm.*m.*g);
    
    % (un)normalized lateral weight transfer
    LatLT  = (LatLT_norm.*m.*g);
    
    % Staitc loads at the tires
    Fz10 = (m*g*b)/(2*wheelbase) ; % static load 
    Fz20 = (m*g*b)/(2*wheelbase) ; % static load 
    Fz30 = (m*g*a)/(2*wheelbase) ; % static load   
    Fz40 = (m*g*a)/(2*wheelbase) ; % static load    
    
    % Downforce
    downforce       = 0.5*airDens*liftCoeff*frontalArea*vx.^2;
    frontDownforce  = aeroBalance*downforce;
    rearDownforce   = (1-aeroBalance)*downforce;
    
    % Lateral LT
    deltaFzf_lat   = LatLT.*LatLTD;   
    deltaFzr_lat   = LatLT.*(1-LatLTD);   

    % Longitudinal WT
    deltaFzf_long   = LongLT*0.5;
    deltaFzr_long   = LongLT*0.5;
    
    % Dynamic Tire Loads
    Fz1 = (Fz10 -deltaFzf_long -deltaFzf_lat + frontDownforce/2 );
    Fz2 = (Fz20 -deltaFzf_long +deltaFzf_lat + frontDownforce/2 );
    Fz3 = (Fz30 +deltaFzr_long -deltaFzr_lat + rearDownforce/2 );
    Fz4 = (Fz40 +deltaFzr_long +deltaFzr_lat + rearDownforce/2 );
       

    %% Tire Slip Calculation
       
    alpha1 = -atan( (cos(delta).*(r.*a+vy) - sin(delta).*(r.*frontTrack./2 + vx)) ./...
                    (cos(delta).*(r.*frontTrack/2 + vx) + sin(delta).*(r.*a + vy) ) ) ;
                
    alpha2 = -atan( (sin(delta).*(r.*frontTrack./2-vx) + cos(delta).*(r.*a + vy)) ./...
                    (cos(delta).*(-r.*frontTrack/2+vx) + sin(delta).*(r.*a+vy)) ) ;
                
    alpha3 = -atan( (vy - r.*b) ./ (vx+r.*rearTrack/2) ) ;
    
    alpha4 = -atan( (vy - r.*b) ./ (vx-r.*rearTrack/2) ) ;
    
    Vxfl = cos(delta).*(r.*frontTrack/2 + vx) + sin(delta).*(r.*a + vy);
    Vyfl = cos(delta).*(r.*a+vy) - sin(delta).*(r.*frontTrack./2 + vx);
    
    Vxfr = cos(delta).*(-r.*frontTrack/2+vx) + sin(delta).*(r.*a+vy);
    Vyfr = sin(delta).*(r.*frontTrack./2-vx) + cos(delta).*(r.*a + vy);
    
    Vxrl = (vx+r.*rearTrack/2);
    Vyrl = (vy - r.*b);
    
    Vxrr = (vx-r.*rearTrack/2);
    Vyrr = (vy - r.*b);
    
    %% slip ratio calculation    

    kappa1 =  -( 1 - (Rl1*omega1./Vxfl) ) ;
    kappa2 =  -( 1 - (Rl2*omega2./Vxfr) ) ;
    kappa3 =  -( 1 - (Rl3*omega3./Vxrl) ) ;
    kappa4 =  -( 1 - (Rl4*omega4./Vxrr) ) ;


    %% external forces at the tire

    [Fx1,Fy1,Mz1]  = Solver.MF5ss_eval(Fz1,kappa1,alpha1,0,MF_f);    
    [Fx2,Fy2,Mz2]  = Solver.MF5ss_eval(Fz2,kappa2,alpha2,0,MF_f);
    [Fx3,Fy3,Mz3]  = Solver.MF5ss_eval(Fz3,kappa3,alpha3,0,MF_r);    
    [Fx4,Fy4,Mz4]  = Solver.MF5ss_eval(Fz4,kappa4,alpha4,0,MF_r);   

      
    % Drag force
    dragForce = 0.5*airDens*dragCoeff*frontalArea*vx.^2;    
    
    % Sum of the forces
    Fx  = (Fx1+Fx2).*cos(delta) - (Fy1+Fy2).*sin(delta) + (Fx3+Fx4) - dragForce ;
    
    Fyf = (Fx1+Fx2).*sin(delta) + (Fy1+Fy2).*cos(delta);
    Fyr = (Fy3+Fy4);
    
    Fy  = Fyf+Fyr ;        
    
    YMfromFy = a*( (Fx1+Fx2).*sin(delta) + (Fy1+Fy2).*cos(delta) ) - b*(Fy3+Fy4);  
    YMfromFx = frontTrack/2.*((Fx1-Fx2).*cos(delta) - (Fy1-Fy2).*sin(delta)) + rearTrack/2.*(Fx3-Fx4);
    YMfromMz = Mz1+Mz2+Mz3+Mz4;

    YMnet = YMfromFy+YMfromFx+YMfromMz;
    
    %% Wheel torques

%     [Tfl,Tfr,Trl,Trr,~,~] = Solver.get_wheel_torque3(tau,Vehicle,x,p);
    [Tfl,Tfr,Trl,Trr,~,~] = Solver.getWheelTorque2(tau,p,Vehicle,x);
    
    %% States Time Derivative
    
    dx = [ (Fx./m) + (r.*vy);
     (Fy./m) - (r.*vx);
     (YMfromFy + YMfromFx + YMfromMz) ./ Iz;
     (Tfl - (Rl1.*Fx1) )./ Iwheel_1;
     (Tfr - (Fx2.*Rl2) )./ Iwheel_2;
     (Trl - (Fx3.*Rl3) )./ Iwheel_3;
     (Trr - (Fx4.*Rl4) )./ Iwheel_4;
     steerRate;
     tauRate ];
       
    %% CONSTRAINTS 
    
    %% Equality Constraints : Load Transfer

    q1  = ((Fx.*CoG)./(wheelbase.*m.*g)) - LongLT_norm;
    q2  = ((Fy.*CoG)./(0.5*(frontTrack+rearTrack).*m.*g)) - LatLT_norm;
    
    % Group them into a matrix    
    q_eq = [q1;
           q2];

    %% Equality Constraints : Tire Saturation
    epsKap = 1e-5;
    epsAlp = 1e-5;
    
    [Fx1_t,~,~]     = Solver.MF5ss_eval(Fz1,kappa1 + epsKap,alpha1,0,MF_f);
    [Fx2_t,~,~]     = Solver.MF5ss_eval(Fz2,kappa2 + epsKap,alpha2,0,MF_f);    
    [Fx3_t,~,~]     = Solver.MF5ss_eval(Fz3,kappa3 + epsKap,alpha3,0,MF_r);
    [Fx4_t,~,~]     = Solver.MF5ss_eval(Fz4,kappa4 + epsKap,alpha4,0,MF_r);   
    
    [~,Fy1_t,~]     = Solver.MF5ss_eval(Fz1, kappa1, alpha1 + epsAlp,0, MF_f);
    [~,Fy2_t,~]     = Solver.MF5ss_eval(Fz2, kappa2, alpha2 + epsAlp,0, MF_f);
    [~,Fy3_t,~]     = Solver.MF5ss_eval(Fz3, kappa3, alpha3 + epsAlp,0, MF_r);
    [~,Fy4_t,~]     = Solver.MF5ss_eval(Fz4, kappa4, alpha4 + epsAlp,0, MF_r);

    
    % Calculation of Slip Stiffness
    C_x_1 = (Fx1_t - Fx1)/epsKap;
    C_x_2 = (Fx2_t - Fx2)/epsKap;
    C_x_3 = (Fx3_t - Fx3)/epsKap;
    C_x_4 = (Fx4_t - Fx4)/epsKap;
    
    % Calculation of Cornering Stiffness
    C_y_1 = (Fy1_t - Fy1)/epsAlp;
    C_y_2 = (Fy2_t - Fy2)/epsAlp;
    C_y_3 = (Fy3_t - Fy3)/epsAlp;
    C_y_4 = (Fy4_t - Fy4)/epsAlp;   

    % Groups the inequality constraint into a single vector
%     g_ineq = [C_x_1;
%      C_x_2;
%      C_y_1;
%      C_y_2;
%      C_y_3;
%      C_y_4;
%      C_x_3;
%      C_x_4];

    %% Equality Constraints : Vehicle Stability

    % Understeer angle
    wb = a+b;
    delta_kin = wb*r./vx;
    theta_uang = delta-delta_kin;
    theta_lim_uang = deg2rad(8);
    C_uang = (theta_uang/theta_lim_uang).^2-1;

    % Bounding the maximum value of beta
    beta = atan(vy./vx);
    beta_lim = deg2rad(5);
    C_beta = (beta/beta_lim).^2-1;

%     yaw_stiffness = jacobian(YMnet,beta);   %yaw stiffness
%     gradient(dot(A,A),A)
%         yaw_stiffness = evalf(gradient(YMnet))./evalf(gradient(beta));   %yaw stiffness
    
%     g_ineq = -[C_x_1;
%      C_x_2;
%      C_x_3;
%      C_x_4];

    g_ineq = -[C_x_1;
     C_x_2;
     C_x_3;
     C_x_4;...
     -C_uang;...
     -C_beta];

 

    %% Saving Constraints : Fuel Usage

    mass_flow  = Solver.mass_flow_calculation(x,tau,Vehicle);
    fuel_usage = ds*(Sf.*mass_flow)*ones(length(n),1) ;

    %% Saving Constraints : Tire Energy
      
    % Longitudinal Sliding Speed
    long_sliding_speed_1 = Vxfl.*kappa1;
    long_sliding_speed_2 = Vxfr.*kappa2;
    long_sliding_speed_3 = Vxrl.*kappa3;
    long_sliding_speed_4 = Vxrr.*kappa4;

    % Lateral Sliding Speed
    lat_sliding_speed_1 = Vxfl.*tan(alpha1);
    lat_sliding_speed_2 = Vxfr.*tan(alpha2);
    lat_sliding_speed_3 = Vxrl.*tan(alpha3);
    lat_sliding_speed_4 = Vxrr.*tan(alpha4);

    % IMPORTANT
    % The way the tire energies are being calculated are not entirely
    % correct now, because I'm taking the integral over the tire power over
    % the distance instead of time. Because the speed is one of the states
    % of the problem and therefore it is a casadi MX variable, I cannot
    % calculate the (variable) time step

    sliding_power_lateral_1 = lat_sliding_speed_1 .* Fy1 ; %[w]
%     sliding_energy_lateral_1 = cumtrapz(cumLaptime,abs(sliding_power_lateral_1))./1; % J
    sliding_energy_lateral_1 = ds*abs(sliding_power_lateral_1)*ones(length(n),1); % J

    sliding_power_lateral_2 = lat_sliding_speed_2 .* Fy2 ; %[w]
%     sliding_energy_lateral_2 = cumtrapz(cumLaptime,abs(sliding_power_lateral_2))./1; % J
    sliding_energy_lateral_2 = ds*abs(sliding_power_lateral_2)*ones(length(n),1); % J

    sliding_power_lateral_3 = lat_sliding_speed_3 .* Fy3 ; %[w]
%     sliding_energy_lateral_3 = cumtrapz(cumLaptime,abs(sliding_power_lateral_3))./1; % J
    sliding_energy_lateral_3 = ds*abs(sliding_power_lateral_3)*ones(length(n),1); % J

    sliding_power_lateral_4 = lat_sliding_speed_4 .* Fy4 ; %[w]
%     sliding_energy_lateral_4 = cumtrapz(cumLaptime,abs(sliding_power_lateral_4))./1; % J 
    sliding_energy_lateral_4 = ds*abs(sliding_power_lateral_4)*ones(length(n),1); % J 

    sliding_power_longitudinal_1 = long_sliding_speed_1 .* Fx1 ; %[w]
    sliding_energy_longitudinal_1 = ds*abs(sliding_power_longitudinal_1)*ones(length(n),1); % J

    sliding_power_longitudinal_2 = long_sliding_speed_2 .* Fx2 ; %[w]
    sliding_energy_longitudinal_2 = ds*abs(sliding_power_longitudinal_2)*ones(length(n),1); % J

    sliding_power_longitudinal_3 = long_sliding_speed_3 .* Fx3 ; %[w]
    sliding_energy_longitudinal_3 = ds*abs(sliding_power_longitudinal_3)*ones(length(n),1); % J
    
    sliding_power_longitudinal_4 = long_sliding_speed_4 .* Fx4 ; %[w]
    sliding_energy_longitudinal_4 = ds*abs(sliding_power_longitudinal_4)*ones(length(n),1); % J

    % Combined Energy Dissipated into the 4 Tires

    fl_tire_energy_combined = sliding_energy_lateral_1 + ...
        sliding_energy_longitudinal_1;

    fr_tire_energy_combined = sliding_energy_lateral_2 + ...
        sliding_energy_longitudinal_2;

    rl_tire_energy_combined = sliding_energy_lateral_3 + ...
        sliding_energy_longitudinal_3;

    rr_tire_energy_combined = sliding_energy_lateral_4 + ...
        sliding_energy_longitudinal_4;    
 
  

    %% saving constraints  


%     saving_constraints.tire_energy = fl_tire_energy_combined + fr_tire_energy_combined +...
%         rl_tire_energy_combined + rr_tire_energy_combined;

    saving_constraints.tire_energy = fl_tire_energy_combined;

    saving_constraints.fuel = fuel_usage(end);

    


end