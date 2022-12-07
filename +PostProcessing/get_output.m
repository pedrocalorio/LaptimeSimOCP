function O = get_output(x,u,p,Vehicle,Track)

%     nGrid = length(Track.sLap);

    kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');
    % Step
    ds = (Track.sLap(end) - Track.sLap(1)) / (length(x(1,:)) - 1);
    
    % Gravity constant
    g = 9.81;
    
    %% System states
    
    n       = x(1,:);
    xi      = x(2,:);
    vx      = x(3,:); % velocity
    vy      = x(4,:); % body-slip angle
    r       = x(5,:); %yaw rate
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

    Sf = (1 - n.*kappa)./(vx.*cos(xi)-vy.*sin(xi));
%     cumLaptime = ds*Sf*ones(length(n),1);
    cumLaptime =cumtrapz(Track.sLap,Sf);
    
    %% Design parameters of the car
    
    %aero
    liftCoeff   = Vehicle.aero.liftCoeff;
    dragCoeff   = Vehicle.aero.dragCoeff;
    frontalArea = Vehicle.aero.frontalArea;
    airDens     = Vehicle.aero.airDensity;
    aeroBalance = p(2);
    
    % susp
    LatLTD  = p(3);
    
    % steering
    steringRatio = Vehicle.steeringratio;
    
    %chassis
    m           = Vehicle.chassis.Ms;
    wheelbase   = Vehicle.chassis.wheelbase;
    weight_dist = p(1);
    a           = wheelbase * (1-weight_dist);
    b           = wheelbase * weight_dist;
    frontTrack  = Vehicle.chassis.frontTrack;
    rearTrack   = Vehicle.chassis.rearTrack;
    Iz          = Vehicle.chassis.Iz;
    CoG         = Vehicle.chassis.CoG;
    

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
    
    % steering angles at the tires
    delta = (steer/steringRatio) ; %radians  
        
    % unnormalized longitudinal weight transfer
    LongLT = (LongLT_norm.*m.*g);
    
    % unnormalized lateral weight transfer
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
    
%     % Lateral LT
    deltaFzf_lat   = LatLT.*LatLTD;   
    deltaFzr_lat   = LatLT.*(1-LatLTD);   

%     % Longitudinal WT
    deltaFzf_long   = LongLT*0.5;
    deltaFzr_long   = LongLT*0.5;
    
    % Tire Loads with downforce
    Fz1 = (Fz10 -deltaFzf_long -deltaFzf_lat + frontDownforce/2 );
    Fz2 = (Fz20 -deltaFzf_long +deltaFzf_lat + frontDownforce/2 );
    Fz3 = (Fz30 +deltaFzr_long -deltaFzr_lat + rearDownforce/2 );
    Fz4 = (Fz40 +deltaFzr_long +deltaFzr_lat + rearDownforce/2 );
       

    %% slip angle calculation
       
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

    % WHEN USING MFEVAL TO CALCULATE TIRE FORCES FROM MAGIC FORMULA

%     side_left = 1;
%     side_right = -1;
%     
%     arrayZeros = zeros(1,length(kappa1))';
%     
%     mfInputs_1 = [Fz1',kappa1',-alpha1'*side_left,...
%         arrayZeros,arrayZeros,omega1'.*Rl1,1.8e5*ones(length(kappa1),1),omega1'];
%     mfInputs_2 = [Fz2',kappa2',-alpha2'*side_right,...
%         arrayZeros,arrayZeros,omega2'.*Rl2,1.8e5*ones(length(kappa1),1),omega2'];
%     mfInputs_3 = [Fz3',kappa3',-alpha3'*side_left,...
%         arrayZeros,arrayZeros,omega3'.*Rl3,1.8e5*ones(length(kappa1),1),omega3'];
%     mfInputs_4 = [Fz4',kappa4',-alpha4'*side_right,...
%         arrayZeros,arrayZeros,omega4'.*Rl4,1.8e5*ones(length(kappa1),1),omega4'];
%     
%     mfeval_o_1 = mfeval(Vehicle.tire_1.MF, mfInputs_1, 211);   
%     mfeval_o_2 = mfeval(Vehicle.tire_2.MF, mfInputs_2, 211);   
%     mfeval_o_3 = mfeval(Vehicle.tire_3.MF, mfInputs_3, 211);   
%     mfeval_o_4 = mfeval(Vehicle.tire_4.MF, mfInputs_4, 211);   
%     
%     [Fx1,Fy1,Mz1,C_x_1,C_y_1] = get_tire_cp_forces(mfeval_o_1, side_left);
%     [Fx2,Fy2,Mz2,C_x_2,C_y_2] = get_tire_cp_forces(mfeval_o_2, side_right);
%     [Fx3,Fy3,Mz3,C_x_3,C_y_3] = get_tire_cp_forces(mfeval_o_3, side_left);
%     [Fx4,Fy4,Mz4,C_x_4,C_y_4] = get_tire_cp_forces(mfeval_o_4, side_right);

    
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
    
    %% Wheel torques

%     [Tfl,Tfr,Trl,Trr,throttle,brakes] = Solver.get_wheel_torque3(tau,Vehicle,x,p);
    [Tfl,Tfr,Trl,Trr,throttle,brakes] = Solver.getWheelTorque2(tau,p,Vehicle,x);
    
    %% Calculates the states TIME derivative
    
%     dx(1,:) =  (Fx./m) + (r.*vy);
%     dx(2,:) =  (Fy./m) - (r.*vx);
%     dx(3,:) =  (YMfromFy + YMfromFx + YMfromMz) ./ Iz;
%     dx(4,:) = (Tfl - (Rl1.*Fx1) )./ Iwheel_1;
%     dx(5,:) = (Tfr - (Fx2.*Rl2) )./ Iwheel_2;
%     dx(6,:) = (Trl - (Fx3.*Rl3) )./ Iwheel_3;
%     dx(7,:) = (Trr - (Fx4.*Rl4) )./ Iwheel_4;
%     dx(8,:) = steerRate;
%     dx(9,:) = tauRate;
%     dx = [ (Fx./m) + (r.*vy);
%      (Fy./m) - (r.*vx);
%      (YMfromFy + YMfromFx + YMfromMz) ./ Iz;
%      (Tfl - (Rl1.*Fx1) )./ Iwheel_1;
%      (Tfr - (Fx2.*Rl2) )./ Iwheel_2;
%      (Trl - (Fx3.*Rl3) )./ Iwheel_3;
%      (Trr - (Fx4.*Rl4) )./ Iwheel_4;
%      steerRate;
%      tauRate ];
       
        %% output for the constraints  
    
    % equality constraints on load transfer
%     q1  = ((Fx.*CoG)./(wheelbase.*m.*g)) - LongLT_norm;
%     q2  = ((Fy.*CoG)./(0.5*(frontTrack+rearTrack).*m.*g)) - LatLT_norm;
    
%     
%     q_eq = [q1;
%            q2];

        % Inequality constraints on tire saturation
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
    
    % g2 and g3
    C_x_1 = (Fx1_t - Fx1)/epsKap;
    C_x_2 = (Fx2_t - Fx2)/epsKap;
    C_x_3 = (Fx3_t - Fx3)/epsKap;
    C_x_4 = (Fx4_t - Fx4)/epsKap;
    
    % additional constraints for lateral tire sat
    C_y_1 = (Fy1_t - Fy1)/epsAlp;
    C_y_2 = (Fy2_t - Fy2)/epsAlp;
    C_y_3 = (Fy3_t - Fy3)/epsAlp;
    C_y_4 = (Fy4_t - Fy4)/epsAlp;   

    % Inq constraints for stability index
    Cf = (C_y_1 + C_y_2);
    Cr = (C_y_3 + C_y_4);

    A_mat = zeros(2,2,length(vx));
    eigenvalues = zeros(2,length(vx));

    for ii=1:length(vx)

        A_mat(:,:,ii) = - [(Cf(ii)+Cr(ii))./(m*vx(ii)) vx(ii)+((a*Cf(ii)-b*Cr(ii))./(m*vx(ii)));...
             (a*Cf(ii)-b*Cr(ii))./(Iz*vx(ii)) (a^2*Cf(ii) + b^2*Cr(ii))./(Iz*vx(ii))];

        eigenvalues(:,ii) = real(eig(A_mat(:,:,ii)));

    end

 



    

    % Inequality constraint on fuel saving 
    we = (omega3 + omega4)/2;
    fTauPos = tau .* (tanh(100*tau) + 1)/2;
    Te = fTauPos.*( Vehicle.engine.maximum_power./we ); 
    n_diesel = 0.40; % combustion efficiency
    e_diesel = 45e6; % energy density [J/kg] 
    
    mass_flow = (1/(n_diesel*e_diesel)).*Te.*we; %[kg/s]
    
%     fuel_usage = trapz(Track.sLap, (Sf.*mass_flow) ) ;
    fuel_usage_int = ds*(Sf.*mass_flow)*ones(length(n),1) ;
    fuel_usage = fuel_usage_int(end);

      
    long_sliding_speed_1 = Vxfl.*kappa1;
    long_sliding_speed_2 = Vxfr.*kappa2;
    long_sliding_speed_3 = Vxrl.*kappa3;
    long_sliding_speed_4 = Vxrr.*kappa4;

    lat_sliding_speed_1 = Vxfl.*tan(alpha1);
    lat_sliding_speed_2 = Vxfr.*tan(alpha2);
    lat_sliding_speed_3 = Vxrl.*tan(alpha3);
    lat_sliding_speed_4 = Vxrr.*tan(alpha4);

    sliding_power_lateral_1 = lat_sliding_speed_1 .* Fy1 ; %[w]
    sliding_energy_lateral_1 = cumtrapz(cumLaptime,abs(sliding_power_lateral_1))./1; % J
%     sliding_energy_lateral_1 = ds*abs(sliding_power_lateral_1)*ones(length(n),1); % J

    sliding_power_lateral_2 = lat_sliding_speed_2 .* Fy2 ; %[w]
    sliding_energy_lateral_2 = cumtrapz(cumLaptime,abs(sliding_power_lateral_2))./1; % J
%     sliding_energy_lateral_2 = ds*abs(sliding_power_lateral_2)*ones(length(n),1); % J

    sliding_power_lateral_3 = lat_sliding_speed_3 .* Fy3 ; %[w]
    sliding_energy_lateral_3 = cumtrapz(cumLaptime,abs(sliding_power_lateral_3))./1; % J
%     sliding_energy_lateral_3 = ds*abs(sliding_power_lateral_3)*ones(length(n),1); % J

    sliding_power_lateral_4 = lat_sliding_speed_4 .* Fy4 ; %[w]
    sliding_energy_lateral_4 = cumtrapz(cumLaptime,abs(sliding_power_lateral_4))./1; % J 
%     sliding_energy_lateral_4 = ds*abs(sliding_power_lateral_4)*ones(length(n),1); % J 

    sliding_power_longitudinal_1 = long_sliding_speed_1 .* Fx1 ; %[w]
    sliding_energy_longitudinal_1 = cumtrapz(cumLaptime,abs(sliding_power_longitudinal_1))./1; % J

    sliding_power_longitudinal_2 = long_sliding_speed_2 .* Fx2 ; %[w]
    sliding_energy_longitudinal_2 = cumtrapz(cumLaptime,abs(sliding_power_longitudinal_2))./1; % J

    sliding_power_longitudinal_3 = long_sliding_speed_3 .* Fx3 ; %[w]
    sliding_energy_longitudinal_3 = cumtrapz(cumLaptime,abs(sliding_power_longitudinal_3))./1; % J
    
    sliding_power_longitudinal_4 = long_sliding_speed_4 .* Fx4 ; %[w]
    sliding_energy_longitudinal_4 = cumtrapz(cumLaptime,abs(sliding_power_longitudinal_4))./1; % J
 
%     O = [];
    thetaModel = interp1(Track.distance,Track.aYaw,Track.sLap) + x(2,:);
    
    % xModel = interp1(Track.distance,Track.xCar,Track.sLap) - x_star(1,:).*sin(thetaModel);
    % yModel = interp1(Track.distance,Track.yCar,Track.sLap) + x_star(1,:).*cos(thetaModel);
    
    xModel = interp1(Track.distance,Track.xCar,Track.sLap) + x(1,:).*sin(thetaModel);
    yModel = interp1(Track.distance,Track.yCar,Track.sLap) - x(1,:).*cos(thetaModel);
    
    O = zeros(56,length(x(1,:))); 
    O(1,:) = cumLaptime;
    O(2,:) = (Fy./m) ;
    O(3,:) = deltaFzf_lat;
    O(4,:) = deltaFzr_lat;
    O(5,:) = deltaFzf_long;
    O(6,:) = deltaFzr_long;
    O(7,:) = throttle;
    O(8,:) = brakes;
    O(9,:) = vx;
    O(10,:) = sliding_energy_lateral_1;
    O(11,:) = sliding_energy_lateral_2;
    O(12,:) = sliding_energy_lateral_3;
    O(13,:) = sliding_energy_lateral_4;
    O(14,:) = sliding_energy_longitudinal_1;
    O(15,:) = sliding_energy_longitudinal_2;
    O(16,:) = sliding_energy_longitudinal_3;
    O(17,:) = sliding_energy_longitudinal_4;
    O(18,:) = Fz1;
    O(19,:) = Fz2;
    O(20,:) = Fz3;
    O(21,:) = Fz4;
    O(22,:) = Fy1;
    O(23,:) = Fy2;
    O(24,:) = Fy3;
    O(25,:) = Fy4;
    O(26,:) = Fx1;
    O(27,:) = Fx2;
    O(28,:) = Fx3;
    O(29,:) = Fx4;
    O(30,:) = alpha1;
    O(31,:) = alpha2;
    O(32,:) = alpha3;
    O(33,:) = alpha4;
    O(34,:) = kappa1;
    O(35,:) = kappa2;
    O(36,:) = kappa3;
    O(37,:) = kappa4;
    O(38,:) = YMfromFx;
    O(39,:) = YMfromFy;
    O(40,:) = YMfromMz;
    O(41,:) = C_x_1;
    O(42,:) = C_y_1;
    O(43,:) = C_y_2;
    O(44,:) = mass_flow;
    O(45,:) = Fx/m;
    O(46,:) = sliding_power_lateral_1;
    O(47,:) = sliding_power_lateral_2;
    O(48,:) = sliding_power_lateral_3;
    O(49,:) = sliding_power_lateral_4;
    O(50,:) = sliding_power_longitudinal_1;
    O(51,:) = sliding_power_longitudinal_2;
    O(52,:) = sliding_power_longitudinal_3;
    O(53,:) = sliding_power_longitudinal_4;
    O(54,:) = rad2deg(steer);
    O(55,:) = xModel;
    O(56,:) = yModel;
 

end