function f = fnObjective(x,Track)
% 
%  f = f(x,y) = integrand of the objective function
%

% vehicle states
n       = x(1,:); % n - distance perpendicular to center line
zeta    = x(2,:); % zeta - angle relative to center line
vx      = x(3,:); % velocity
vy      = x(4,:); % body-slip angle
% r       = x(5,:); %yaw rate

kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

% conversion from time to distance
Sf = (1 - n.*kappa)./(vx.*cos(zeta)-vy.*sin(zeta));

f = Sf;

end