function [fWheelTorque, rWheelTorque] = getWheelTorque( u, x, vehicle )

%     fTau = u;
%     nWheel = x;
% 
%     PWR_MaxPower = vehicle.engine.maximum_power;
%     PWR_MaxBrkPower = 0;
%     PWR_MaxBrkTorque = vehicle.brakes.maximum_torque;
%   
%     %{
%     if fTau > 0
%         fTractiveTorqueF = 0.0;
%         fTractiveTorqueR = PWR_MaxPower * fTau / nWheel;
%     else
%         fTractiveTorqueF = PWR_MaxBrkTorque * fTau;
%         fTractiveTorqueR = PWR_MaxBrkTorque * fTau;
%     end
%     %}
%     
%     % The code below is replacing the if statements above, but smoothing the
%     % transition between throttle and braking.
%     fTauPos = fTau .* (tanh(20*fTau) + 1)/2;
%     fTauNeg = fTau .* (tanh(-20*fTau) + 1)/2;
%     
%     fWheelTorque = ( (PWR_MaxPower + PWR_MaxBrkPower) * fTauPos - PWR_MaxBrkPower ) ./ nWheel;
%     fWheelTorque = fWheelTorque + PWR_MaxBrkTorque * fTauNeg;
    
    



% braking_torque = u * vehicle.brakes.maximum_torque;
% Torque_front = braking_torque*vehicle.brakes.bias;
% Torque_rear  = braking_torque*(1-vehicle.brakes.bias);

Torque_front = zeros(1,length(u));
Torque_rear  = zeros(1,length(u));


fTauPos = u .* (tanh(20*u) + 1)/2;
fTauNeg = u .* (tanh(-20*u) + 1)/2;

for i=1:length(u)

    if u(i) > 0
        engine = vehicle.engine;

        % no wheel torque in the fronts when the car is acceleratig
        Torque_front(i) = 0;

        % calculates engine rotational velocity    
        engineRotVel = ( x(4,i)+x(5,i) )/2 * engine.Gr*30/pi;

        Te = engine.maximum_power * interp1(engine.RPM, engine.Power_R, engineRotVel./9.549,'pchip','extrap') * ...
            interp1(engine.Throttle, engine.Power_T, 100*u(i)) ;

        Torque_rear(i) = Te*engine.Gr/2; 
    else
        %
        braking_torque = u(i) * vehicle.brakes.maximum_torque;
        Torque_front(i) = braking_torque*vehicle.brakes.bias/2;
        Torque_rear(i)  = braking_torque*(1-vehicle.brakes.bias)/2;
    end

end

fWheelTorque=Torque_front;
rWheelTorque=Torque_rear;

end