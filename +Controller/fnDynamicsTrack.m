function dx = fnDynamicsTrack(x,Track)
addpath('C:/dev/libraries/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

n       = x(1,:); % n - distance perpendicular to center line
zeta    = x(2,:); % zeta - angle relative to center line
vx      = x(3,:); % velocity
vy      = x(4,:); % body-slip angle
r       = x(5,:); %yaw rate

Sf = (1 - n.*kappa)./(vx.*cos(zeta)-vy.*sin(zeta));

dx = [Sf .* ( vx.*sin(zeta) + vy.*cos(zeta) );
        (Sf .* r) - kappa];