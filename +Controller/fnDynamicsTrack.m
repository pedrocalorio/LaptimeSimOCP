function dx = fnDynamicsTrack(x,states,Track)

kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

n       = x(1,:); % velocity
zeta    = x(2,:); % body-slip angle

vx      = states(1,:); % velocity
vy      = states(2,:); % body-slip angle
r       = states(3,:); % body-slip angle

Sf = (1 - n.*kappa)./(vx.*cos(zeta)-vy.*sin(zeta));

dx = zeros(2,length(n));

dx(1,:) = Sf .* ( vx.*sin(zeta) + vy.*cos(zeta) );
dx(2,:) = (Sf .* r) - kappa;