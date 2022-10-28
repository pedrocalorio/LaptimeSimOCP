function soln = fnDirectCollocation(problem)

%% To make code more readable

B       = problem.bounds;
F       = problem.func;
Opt     = problem.options;
sLap    = problem.dsSystem.td.sLap;
nGrid   = length(F.weights);
nDesign = length(problem.bounds.design.low);
Track   = problem.dsSystem.td;
Vehicle = problem.dsSystem.vd;

savingConstraintsBounds.fuel = problem.bounds.fuel_consumption;
savingConstraintsBounds.tire_energy = problem.bounds.tire_energy;

%%  Pack the initial guess

[zGuess, pack] = Utilities.packDecVar_mex(problem.guess.state, problem.guess.control);

if nDesign~=0
    zGuess = [zGuess; problem.guess.design];
end


%%  Pack the lower and upper bounds

xLow = [B.initialState.low, B.state.low*ones(1,nGrid-2), B.finalState.low];
uLow = [B.initialControl.low, B.control.low*ones(1,nGrid-2), B.finalControl.low];
% Special condition for when the track is generated with GPS
xLow(1,:) = interp1(Track.distance,Track.left_offset,Track.sLap) + Vehicle.chassis.frontTrack;
% Joins the state and control variables into one
zLow = Utilities.packDecVar_mex(xLow,uLow);

if nDesign~=0
    zLow = [zLow; B.design.low];
end

xUpp = [B.initialState.upp, B.state.upp*ones(1,nGrid-2), B.finalState.upp];
uUpp = [B.initialControl.upp, B.control.upp*ones(1,nGrid-2), B.finalControl.upp];
% Special condition for when the track is generated with GPS
xUpp(1,:) = interp1(Track.distance,Track.right_offset,Track.sLap) - Vehicle.chassis.frontTrack;
% Joins the state and control variables into one
zUpp = Utilities.packDecVar_mex(xUpp,uUpp);

if nDesign~=0
    zUpp = [zUpp; B.design.upp];
end

%% Set up the functions, bounds, and options for fmincon

P.objective = @(z)( myObjective(z, pack, F.weights, Track, nDesign) );
P.nonlcon   = @(z)( myConstraint(z, pack, F.defectCst, Track, Vehicle, nDesign, savingConstraintsBounds) );

P.x0 = zGuess;
P.lb = zLow;
P.ub = zUpp;
P.Aineq = []; P.bineq = []; % Unused
P.Aeq = []; P.beq = [];     % Unused
P.options = Opt.nlpOpt;
P.solver = 'fmincon';

%% Call fmincon to solve the non-linear program (NLP)
tic;
[zSoln, objVal, exitFlag, output] = fmincon(P);
nlpTime = toc;

%% Unpack the solution and store the results

[xSoln,uSoln] = Utilities.unPackDecVar_mex(zSoln(1:end-nDesign),pack); % Unpack decision variables

soln.grid.sLap = sLap;
soln.grid.state = xSoln;
soln.grid.control = uSoln;

if nDesign~=0
    soln.grid.design = zSoln(end-nDesign+1:end);
else
    soln.grid.design = [Vehicle.chassis.frontMassDist,...
        Vehicle.aero.aeroBalance,...
        Vehicle.susp.LatLTD,...
        Vehicle.brakes.bias]';
end

soln.info = output;
soln.info.nlpTime = nlpTime;
soln.info.exitFlag = exitFlag;

% Get the final laptime
% Step
ds = (sLap(end) - sLap(1)) / (pack.nGrid - 1);

% Inverse of the time derivative of the positon
integrand = Controller.fnObjective(xSoln,Track);   % Calculate the integrand of the cost function

% Calculate laptime via integration
% laptime = ds*integrand*ones(length(integrand),1);  % Integration
laptime = trapz(Track.sLap,integrand);  % Integration

soln.info.laptime = laptime;
soln.info.objVal  = sqrt(objVal);

end

%% Utility Functions


%% Objective function

function cost = myObjective(z,pack,weights,Track,nDesign)
% This function returns the final weighted cost

sLap = Track.sLap;

[x,~] = Utilities.unPackDecVar_mex(z(1:end-nDesign),pack);

% Step
ds = (sLap(end) - sLap(1)) / (pack.nGrid - 1);

% Inverse of the time derivative of the positon
integrand = Controller.fnObjective(x,Track);   % Calculate the integrand of the cost function

% Calculate laptime via integration
laptime = ds*integrand*weights;  % Integration

% absTau = abs(x(11,:));
% sumInvTau = trapz(1./absTau)^2;
% 
cost = laptime ;

% cost = laptime;

end

%% Constraint Function

function [c, ceq] = myConstraint(z,pack,defectCst,Track,Vehicle,nDesign,savingConstraintsBound)
% This function computes the defects along the path
% and then evaluates the user-defined constraint functions

sLap = Track.sLap;

kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

[x,u] = Utilities.unPackDecVar_mex(z(1:end-nDesign),pack);

% Calculate defects along the path
ds = (sLap(end) - sLap(1)) / (pack.nGrid - 1);

% separate the vehicle and track states
vehicle_states = x(3:11,:);
track_states   = x(1:2,:);

% gets the vehicle design variables
if nDesign==0
    vehicle_design_variables = [Vehicle.chassis.frontMassDist,...
        Vehicle.aero.aeroBalance,...
        Vehicle.susp.LatLTD,...
        Vehicle.brakes.bias]';
else
    vehicle_design_variables = z(end-nDesign+1:end);
end


% gets the vehicle states derivative IN TIME

[dx_vehicle,g,q,~,saving_constraints] = Controller.fnDynamicsVehicle([],...
                                            vehicle_states,...
                                            u,...
                                            Vehicle,...
                                            vehicle_design_variables, ...
                                            Track,track_states);


% Calculates the conversion factor for the transformation, which is the
% inverse speed along the center line
Sf = (1 - x(1,:).*kappa)./(x(3,:).*cos(x(2,:))-x(4,:).*sin(x(2,:)));

% Transform vehicle state into space/distance dependent 
dx_vehicle = Sf.*dx_vehicle;

% Gets the track states derivative IN SPACE/DISTANCE
dx_track = Controller.fnDynamicsTrack(track_states,vehicle_states,Track);

% Joins the state vector into one
dx = [dx_track; dx_vehicle];

% Gets the defects to enforce continous dynamics
defects = defectCst(ds,x,dx);

% Call user-defined constraints and combine with defects
[c, ceq] = Solver.fnCollectConstraints(defects,g,q,x,u,...
                                        saving_constraints,savingConstraintsBound);


end
