function [problem] = fnInitBounds(problem)

TrackWidth = problem.dsSystem.td.TrackWidth;

Vehicle = problem.dsSystem.vd;

rl_f = Vehicle.tire_1.radius;
rl_r = Vehicle.tire_3.radius;

frontTrack = Vehicle.chassis.frontTrack;
rearTrack = Vehicle.chassis.rearTrack;

VehicleTrack = (frontTrack+rearTrack)/2;

Vx0 = 350/3.6;

% Upper and lower bound on state
problem.bounds.state.low        = [-TrackWidth/2+VehicleTrack/2 -pi -inf*ones(1,14)  0 -inf.*ones(1,5) -inf.*ones(1,4) zeros(1,4) -pi -1 zeros(1,8)]';
problem.bounds.state.upp        = [ TrackWidth/2-VehicleTrack/2  pi +inf*ones(1,14) inf inf.*ones(1,5) inf.*ones(1,4) inf.*ones(1,4) pi 1 inf*ones(1,8)]';

% Upper and lower bound on the initial state
problem.bounds.initialState.low = [-TrackWidth/2+VehicleTrack/2 -pi -inf*ones(1,14)  0 -inf.*ones(1,5) -inf.*ones(1,4) zeros(1,4) -0.05 0 zeros(1,8)]';
problem.bounds.initialState.upp = [ TrackWidth/2-VehicleTrack/2  pi +inf*ones(1,14) inf inf.*ones(1,5) inf.*ones(1,4) inf.*ones(1,4) 0.05 1 inf*ones(1,8)]';

% Upper and lower bound on the final state
problem.bounds.finalState.low   = problem.bounds.initialState.low;
problem.bounds.finalState.upp   = problem.bounds.initialState.upp;

% Upper and lower bound on the control inputs
problem.bounds.control.low = [-pi/4 -5]';
problem.bounds.control.upp = [ pi/4  5]';

problem.bounds.initialControl.low = [-0.001 -5]';
problem.bounds.initialControl.upp = [ 0.001  5]';

problem.bounds.finalControl.low = [-0.001 -5]';
problem.bounds.finalControl.upp = [ 0.001  5]';

%% defines the lower and upper bounds on the values admissible by the optimizer when optimizing the vehicle parameters 

if problem.options.optimalDesignFlag == true                            

    % [weight distribution, downforce distribution, lateral load transfer distribution, brake torque distribution]
    problem.bounds.design.low = [0.40 0.35 3e3 0.55]';
    problem.bounds.design.upp = [0.50 0.50 13e3 0.80]';
    
else

    problem.bounds.design.low = []';
    problem.bounds.design.upp = []';

end

%% defines the value of the maximum fuel consumption admissible 

if problem.options.fuelSavingFlag == true
    problem.bounds.fuel_consumption = 1.2405*0.75;
else
    problem.bounds.fuel_consumption = -1;
end

%% defines the value of the maximum total tire energy dissipated in a lap admissible 

if problem.options.tireEnergySavingFlag == true
    problem.bounds.tire_energy = 4.996e+07*0.6;
else
    problem.bounds.tire_energy = -1;
end


end