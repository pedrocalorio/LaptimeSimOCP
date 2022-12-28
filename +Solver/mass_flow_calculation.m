function mass_flow = mass_flow_calculation(omega3,omega4,tau,Vehicle)
    % Inequality constraint on fuel saving 

    % Calculates the engine speed
    we      = (omega3 + omega4)/2;
    
    fTauPos = tau .* (tanh(100*tau) + 1)/2;
    Te      = fTauPos.*( Vehicle.engine.maximum_power./we ); 
    n_diesel = 0.40; % combustion efficiency
    e_diesel = 45e6; % energy density [J/kg] 
    
    mass_flow = (1/(n_diesel*e_diesel)).*Te.*we; %[kg/s]
    

end