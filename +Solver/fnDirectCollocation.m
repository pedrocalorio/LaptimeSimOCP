function soln = fnDirectCollocation(problem)

addpath('C:/dev/libraries/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*


% To make code more readable
B       = problem.bounds;
F       = problem.func;
nGrid   = length(F.weights);
nDesign = length(problem.bounds.design.low);
Track   = problem.dsSystem.td;
Vehicle = problem.dsSystem.vd;

savingConstraintsBounds.fuel        = problem.bounds.fuel_consumption;
savingConstraintsBounds.tire_energy = problem.bounds.tire_energy;


%%  Pack the initial guess

xGuess = problem.guess.state;
uGuess = problem.guess.control;

if nDesign~=0
    pGuess = problem.guess.design;
end


%%  Pack the lower and upper bounds

xLow = [B.initialState.low, B.state.low*ones(1,nGrid-2), B.finalState.low];
uLow = [B.initialControl.low, B.control.low*ones(1,nGrid-2), B.finalControl.low];
pLow = B.design.low;

% Special condition for when the track is generated with GPS
% xLow(1,:) = interp1(Track.distance,Track.left_offset,Track.sLap) + Vehicle.chassis.frontTrack/2;


xUpp = [B.initialState.upp, B.state.upp*ones(1,nGrid-2), B.finalState.upp];
uUpp = [B.initialControl.upp, B.control.upp*ones(1,nGrid-2), B.finalControl.upp];
pUpp = B.design.upp;

% Special condition for when the track is generated with GPS
% xUpp(1,:) = interp1(Track.distance,Track.right_offset,Track.sLap) - Vehicle.chassis.frontTrack/2;


%%Create the design variables of the optimization problem in CasADi MX variable 
opti = casadi.Opti();
if nDesign~=0    

    xk = opti.variable(40,nGrid);
    uk = opti.variable(2,nGrid);
    pk = opti.variable(4);
    
    % Set up the functions, bounds, and options for fmincon
    
    J =  myObjective(xk, F.weights, Track) ;
%     [q,g] =  myConstraint(xk, uk, pk, F.defectCst, Track, Vehicle, savingConstraintsBounds) ;
    [q,~] =  myConstraint(xk, uk, pk, F.defectCst, Track, Vehicle, savingConstraintsBounds) ;
%     q.generate('gen.c');
    
    opti.minimize(  J   );
    
    opti.subject_to( q==0 );
%     opti.subject_to( g<=0 );
    opti.subject_to( xLow<=xk<=xUpp );
    opti.subject_to( uLow<=uk<=uUpp );
    opti.subject_to( pLow<=pk<=pUpp );
    opti.set_initial(xk, xGuess);
    opti.set_initial(uk, uGuess);
    opti.set_initial(pk, pGuess);
    
    p_opts = struct('expand',true);
    s_opts = struct('max_iter',1e6);
    opti.solver('ipopt',p_opts,...
                    s_opts);
    
    sol = opti.solve();
else

    xk = opti.variable(11,nGrid);
    uk = opti.variable(4,nGrid);
    pk = [Vehicle.chassis.frontMassDist Vehicle.aero.aeroBalance,...
        Vehicle.susp.LatLTD Vehicle.brakes.bias];
    
    %% Set up the functions, bounds, and options for fmincon
    
    J =  myObjective(xk, uk, F.weights, Track) ;
    [q,g] =  myConstraint(xk, uk, pk, F.defectCst, Track, Vehicle, savingConstraintsBounds) ;
%     g =  myIneqConstraint(xk, uk, pk, F.defectCst, Track, Vehicle, savingConstraintsBounds) ;
    
    opti.minimize(  J   );
    
    opti.subject_to( q==0 );
    opti.subject_to( g<=0 );
    opti.subject_to( xLow<=xk<=xUpp );
    opti.subject_to( uLow<=uk<=uUpp );
    opti.set_initial(xk, xGuess);
    opti.set_initial(uk, uGuess);
    
    p_opts = struct('expand',true);
    s_opts = struct('max_iter',1e6);
    opti.solver('ipopt',p_opts,...
                    s_opts);
    
    sol = opti.solve();

end


soln = struct();
soln.state = sol.value(xk);
soln.control = sol.value(uk);
if nDesign~=0   
    soln.design = sol.value(pk);
else
    soln.design = [Vehicle.chassis.frontMassDist Vehicle.aero.aeroBalance,...
        Vehicle.susp.LatLTD Vehicle.brakes.bias];
end
soln.obj_func = (myObjective(soln.state,soln.control, F.weights, Track));


end




%% Objective function

function cost = myObjective(x,weights,Track)
% This function returns the final weighted cost

sLap = Track.sLap;

nGrid = length(sLap);

% [x,~] = Utilities.unPackDecVar(z(1:end),pack);

% Step
ds = (sLap(end) - sLap(1)) / (nGrid - 1);

% Inverse of the time derivative of the positon
integrand = Controller.fnObjectiveCasadi(x,Track);   % Calculate the integrand of the cost function

% Calculate laptime via integration
laptime = ds*integrand*weights;  % Integration
% laptime = trapz(Track.sLap,integrand);  % Integration

% absTau = abs(x(11,:));
% sumInvTau = trapz(1./absTau)^2;
% 
cost = laptime ;

% cost = laptime;

end

%% Constraint Function

function [ceq,c] = myConstraint(x,u,p,defectCst,Track,Vehicle,savingConstraintsBounds)
% This function computes the defects along the path
% and then evaluates the user-defined constraint functions

sLap = Track.sLap;

nGrid = length(sLap);

% kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');
kappa = interp1(Track.sLap,Track.curv,Track.sLap,'spline');

% Calculate defects along the path
ds = (sLap(end) - sLap(1)) / (nGrid - 1);


% gets the vehicle states derivative IN TIME

[dx_vehicle,g,savingConstraints] = Controller.fnDynamics14DOFVehicle(x,u,p,...
                                            Vehicle,...
                                            Track);

n       = x(1,:); % velocity
zeta    = x(2,:); % body-slip angle
vx      = x(17,:); % velocity
vy      = x(18,:); % body-slip angle

% Calculates the conversion factor for the transformation, which is the
% inverse speed along the center line
Sf = (1 - n.*kappa)./(vx.*cos(zeta)-vy.*sin(zeta));

% Transform vehicle state into space/distance dependent 
dx_vehicle = (Sf'.*dx_vehicle')';

% Gets the track states derivative IN SPACE/DISTANCE
dx_track = Controller.fnDynamicsTrack(x,Track);

% Joins the state vector into one
dx = [dx_track; dx_vehicle];

% Gets the defects to enforce continous dynamics
defects = defectCst(ds,x,dx);

% Call user-defined constraints and combine with defects
[ceq,c] = Solver.fnCollectConstraints(defects,g,x,u,savingConstraints,savingConstraintsBounds);



end


