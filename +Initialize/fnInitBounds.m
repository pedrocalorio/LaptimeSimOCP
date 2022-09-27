function [problem] = fnInitBounds(problem)

TrackWidth = problem.dsSystem.td.TrackWidth;

Vehicle = problem.dsSystem.vd;

rl_f = Vehicle.tire_1.radius;
rl_r = Vehicle.tire_3.radius;

frontTrack = Vehicle.chassis.frontTrack;
rearTrack = Vehicle.chassis.rearTrack;

VehicleTrack = (frontTrack+rearTrack)/2;

Vx0 = inf/3.6;

% Upper and lower bound on state
problem.bounds.state.low        = [-TrackWidth/2+VehicleTrack/2 -inf   0 -inf -inf zeros(1,4) -pi -1]';
problem.bounds.state.upp        = [ TrackWidth/2-VehicleTrack/2  inf inf  inf  inf inf.*ones(1,4) pi 1]';

% Upper and lower bound on the initial state
problem.bounds.initialState.low = [-TrackWidth/2+VehicleTrack/2 -pi/4 0/3.6 -0.5 -1e-5 zeros(1,4) -1e-5 0]';
problem.bounds.initialState.upp = [ TrackWidth/2-VehicleTrack/2  pi/4 Vx0    0.5  1e-5 Vx0/rl_f Vx0/rl_f Vx0/rl_r Vx0/rl_r 1e-5 1]';

% Upper and lower bound on the final state
problem.bounds.finalState.low   = [-TrackWidth/2+VehicleTrack/2 -pi/4 0/3.6 -0.5 -1e-5 zeros(1,4) -1e-5 0]';
problem.bounds.finalState.upp   = [ TrackWidth/2-VehicleTrack/2  pi/4 Vx0    0.5  1e-5 Vx0/rl_f Vx0/rl_f Vx0/rl_r Vx0/rl_r 1e-5 1]';

% Upper and lower bound on the control inputs
problem.bounds.control.low = [-pi/4 -3 -inf -inf]';
problem.bounds.control.upp = [ pi/4  3  inf  inf]';

problem.bounds.initialControl.low = [-0.001 -0.1 -inf -inf]';
problem.bounds.initialControl.upp = [ 0.001  0.1  inf  inf]';

problem.bounds.finalControl.low = [-0.001 -0.1 -inf -inf]';
problem.bounds.finalControl.upp = [ 0.001  0.1  inf  inf]';

if problem.options.optimalDesignFlag == true                            

    problem.bounds.design.low = [0.40 0.35 0.40 0.55]';
    problem.bounds.design.upp = [0.50 0.50 0.75 0.80]';
    
else

    problem.bounds.design.low = []';

    problem.bounds.design.upp = []';

end

if problem.options.fuelSavingFlag == true
    problem.bounds.fuel_consumption = 1.2405*0.75;
else
    problem.bounds.fuel_consumption = -1;
end


end