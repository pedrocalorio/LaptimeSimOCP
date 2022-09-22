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

% conversion from time to distance
Sf = (1 - n.*Track.curv)./(vx.*cos(zeta)-vy.*sin(zeta));

f = Sf;

end