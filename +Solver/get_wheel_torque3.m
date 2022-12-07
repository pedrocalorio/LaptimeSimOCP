function [Torque_FL,Torque_FR,Torque_RL,Torque_RR,fTauPos,fTauNeg,omega_e,...
    Pe] = get_wheel_torque3(tau,Vehicle,x,p)

    % Calculated wheel speeds from the system states
    omega3 = x(8,:); %rad/s
    omega4 = x(9,:); %rad/s

    % Calculates positive and negative longitudinal control input variables
    fTauPos = tau .* (tanh(100*tau) + 1)/2;
    fTauNeg = tau .* (tanh(-100*tau) + 1)/2;
    
    % Gear ratio
%     speed = casadi.DM(length(x(3,:)),1)';
    speed = x(3,:);
    Gr = interp1(Vehicle.engine.Gr(:,2),...
        Vehicle.engine.Gr(:,1),...
        speed,... % interpolate speed in km/h
        "linear","extrap");
    
    % Brakes    
%     brk_bias       = Vehicle.brakes.bias;
    brk_bias       = p(4);
    brk_max_torque = Vehicle.brakes.maximum_torque;
    
    % Engine
    % Maximum power available from powertrain
    P_max       = Vehicle.engine.maximum_power;    
    % Maximum power available for brake engine
    Pe_drag     = Vehicle.engine.power_drag;
    
    % Differential
    A = Vehicle.differential.To;
    B = Vehicle.differential.G;
    
    %% Calculations
    
    % Torque from engine brake
    Te_drag = 2*Pe_drag./(omega3+omega4);
    
    % Engine rotational velocity - RPM
    omega_e     = 2.625*Gr.*(omega3+omega4)/2 ; %rad/s
    omega_e_rpm = omega_e * 30/pi ; %rpm
    
    Pe = P_max .* interp1(Vehicle.engine.RPM,...
        Vehicle.engine.Power_R,...
        omega_e_rpm,...
        'spline','extrap') .*...
        interp1(Vehicle.engine.Throttle,...
        Vehicle.engine.Power_T,...
        100.*fTauPos,...
        'spline','extrap');
    
%     Te = fTauPos.*( (2*Pe./(omega3 + omega4)) - Te_drag);   
    Te = fTauPos.*( Pe./omega_e );   
%     Te = fTauPos.*( (Pe./(omega_e)) - Te_drag);   
%     Te = fTauPos.*( (Pe./omega_e) - Te_drag);   
    
    % Calculates the torque difference between the wheels 
    delta_T = A.*sin(atan(B.*(omega4-omega3) ));
    
    Tf     = brk_bias*brk_max_torque .* fTauNeg;
    Tr_brk = ( (1-brk_bias)*brk_max_torque + Te_drag) .* fTauNeg;
    
%     Te =  ( (Pe./ (omega3 + omega4)) - Te_drag) .*fTauPos ;
%     Te =  (Pe * fTauPos ) ./ (omega3 + omega4) ;
    
    % Torque at each wheels
    Torque_FL = Tf/2;
    Torque_FR = Tf/2;
    Torque_RL = ((Te.*2.6250)/2+Tr_brk/2) + delta_T;
    Torque_RR = ((Te.*2.6250)/2+Tr_brk/2) - delta_T;
%     Torque_RL = (Te/2+Tr_brk/2) + delta_T;
%     Torque_RR = (Te/2+Tr_brk/2) - delta_T;

end

