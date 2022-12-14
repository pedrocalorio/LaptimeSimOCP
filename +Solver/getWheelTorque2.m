function [Torque_FL,Torque_FR,Torque_RL,Torque_RR,fTauPos,fTauNeg,...
    omega_e,Te]  = getWheelTorque2(tau,Vehicle,x,p)

    omega3 = x(6,:);
    omega4 = x(7,:);

    % Calculates positive and negative longitudinal control input variables
    fTauPos = tau .* (tanh(100*tau) + 1)/2;
    fTauNeg = tau .* (tanh(-100*tau) + 1)/2;
    
   
    % Brakes    
%     brk_bias       = Vehicle.brakes.bias;
    brk_bias       = p(4);
    brk_max_torque = Vehicle.brakes.maximum_torque;
    
    % Engine
    % maximum power available from powertrain
    P_max       = Vehicle.engine.maximum_power;  
    
    % Differential
    A = Vehicle.differential.To;
    B = Vehicle.differential.G;
    
    %% Calculations 

    Pe = P_max;
    
    omega_e = (omega3 + omega4)/2;
    Te = fTauPos.*( (2*Pe./(omega3 + omega4)) );   
%     Te = fTauPos.*( (Pe./omega_e) - Te_drag);   
%     Te = fTauPos.*( (Pe./(omega_e)) - Te_drag);   
%     Te = fTauPos.*( (Pe./omega_e) - Te_drag);   
    
    delta_T = A.*sin(atan(B.*(x(7,:)-x(6,:)) ));
    
    Tf     = brk_bias*brk_max_torque .* fTauNeg;
    Tr_brk = (1-brk_bias)*brk_max_torque  .* fTauNeg;
    
%     Te =  ( (Pe./ (omega3 + omega4)) - Te_drag) .*fTauPos ;
%     Te =  (Pe * fTauPos ) ./ (omega3 + omega4) ;
    
    % Torque at each wheels
    Torque_FL = Tf/2;
    Torque_FR = Tf/2;
    Torque_RL = (Te/2+Tr_brk/2) + delta_T;
    Torque_RR = (Te/2+Tr_brk/2) - delta_T;


end