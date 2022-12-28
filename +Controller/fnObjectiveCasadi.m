function f = fnObjectiveCasadi(x,Track)
% 
%  f = f(x,y) = integrand of the objective function
%

% vehicle states
n       = x(1,:); % n - distance perpendicular to center line
zeta    = x(2,:); % zeta - angle relative to center line
vx      = x(17,:); % velocity
vy      = x(18,:); % body-slip angle
% r       = x(5,:); %yaw rate

% n       = z(1,(1-1)*nGrid+1:1*nGrid); % n - distance perpendicular to center line
% zeta    = z(1,(2-1)*nGrid+1:2*nGrid); % zeta - angle relative to center line
% vx      = z(1,(3-1)*nGrid+1:3*nGrid); % velocity
% vy      = z(1,(4-1)*nGrid+1:4*nGrid); % body-slip angle

% kappa = interp1(Track.distance,Track.curv,Track.sLap,'linear');
kappa = interp1(Track.sLap,Track.curv,Track.sLap,'linear');

% conversion from time to distance
Sf = (1 - n.*kappa)./(vx.*cos(zeta)-vy.*sin(zeta));

f = Sf;

end