function [Torque_FL,Torque_FR,Torque_RL,Torque_RR,fTauPos,fTauNeg] = get_wheel_torque3(tau,Vehicle,x,p)

    omega3 = x(6,:);
    omega4 = x(7,:);

    % Calculates positive and negative longitudinal control input variables
    fTauPos = tau .* (tanh(100*tau) + 1)/2;
    fTauNeg = tau .* (tanh(-100*tau) + 1)/2;
    
    % Gear ratio
    Gr = Vehicle.engine.Gr;
    
    % Brakes    
%     brk_bias       = Vehicle.brakes.bias;
    brk_bias       = p(4);
    brk_max_torque = Vehicle.brakes.maximum_torque;
    
    % Engine
    % maximum power available from powertrain
    P_max       = Vehicle.engine.maximum_power;    
    % maximum power available for brake engine
    Pe_drag     = Vehicle.engine.power_drag;
    
    % Differential
    A = Vehicle.differential.To;
    B = Vehicle.differential.G;
    
    %% Calculations
    
    % Torque from engine brake
    Te_drag = 2*Pe_drag./(x(6,:)+x(7,:));
    
    % Engine rotational velocity - RPM
    omega_e     = (omega3+omega4)/2 .* Gr;
    omega_e_rpm = omega_e * 30/pi ;
    
    Pe = P_max .* interp1(Vehicle.engine.RPM,Vehicle.engine.Power_R,omega_e_rpm,'spline','extrap') .*...
        interp1(Vehicle.engine.Throttle,Vehicle.engine.Power_T,100.*fTauPos,'spline','extrap');
%     Pe = P_max .* interp1(Vehicle.engine.Throttle,Vehicle.engine.Power_T,100.*fTauPos,'spline','extrap');
    
%     Te = fTauPos.*( (Pe./(omega3 + omega4)) - Te_drag);   
    Te = fTauPos.*( Pe./omega_e );   
%     Te = fTauPos.*( (Pe./(omega_e)) - Te_drag);   
%     Te = fTauPos.*( (Pe./omega_e) - Te_drag);   
    
    delta_T = A.*sin(atan(B.*(x(7,:)-x(6,:)) ));
    
    Tf     = brk_bias*brk_max_torque .* fTauNeg;
    Tr_brk = ( (1-brk_bias)*brk_max_torque + Te_drag) .* fTauNeg;
    
%     Te =  ( (Pe./ (omega3 + omega4)) - Te_drag) .*fTauPos ;
%     Te =  (Pe * fTauPos ) ./ (omega3 + omega4) ;
    
    % Torque at each wheels
    Torque_FL = Tf/2;
    Torque_FR = Tf/2;
    Torque_RL = ((Te*Gr)/2+Tr_brk/2) + delta_T;
    Torque_RR = ((Te*Gr)/2+Tr_brk/2) - delta_T;
%     Torque_RL = (Te/2+Tr_brk/2) + delta_T;
%     Torque_RR = (Te/2+Tr_brk/2) - delta_T;

end

