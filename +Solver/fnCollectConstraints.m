function [c, ceq] = fnCollectConstraints(defects,g,q,x,u,fuel_consumption,fuel_constraint)

% collect deffects from collocation method to put in the ceq matrix of
% fmincon
ceq_dynamics = reshape(defects,numel(defects),1);

% collect the equality constraints on wheel load distribution
ceq_vehicle = reshape(q,numel(q),1);

% initSpeed = x(3,1);
% finalSpeed = x(3,end);

% constraints so that the beggning and the end of the lap is the same
ceq_initial  = [x(1,1) x(2,1) u(3,1) u(4,1)];
% ceq_initial  = [x(1,1) x(2,1) ];
% ceq_initial  = [];
ceq_terminal = [x(1,end) x(2,end) u(3,end) u(4,end)];
% ceq_terminal = [x(1,end) x(2,end) ];
% ceq_terminal = [];

% inequality constraints on having always positve cornering stiffness and
% slip stiffness 
c = reshape(g,numel(g),1);

if fuel_constraint ~= -1
    c = [c; fuel_consumption-fuel_constraint];
end

% inequality constraint so the initial velocity of the lap isn't too far
% from the final 
% c_if = [abs(x(3,1)-x(3,end))-1,... % long vel
%         abs(x(4:5,1)'-x(4:5,end)')-1,... % lat vel and yaw rate
%         abs(x(6:7,1)'-x(6:7,end)')/0.330-1,...% front wheel speeds
%         abs(x(8:9,1)'-x(8:9,end)')/0.330-1,...% rear wheel speeds
%         abs(x(10,1)-x(10,end))-0.05,... % steering
%         abs(x(11,1)-x(11,end))-0.05,... % tau
%         abs(u(3:4,1)'-u(3:4,end)')-0.05]'; % load trans var

c_if = [abs(x(3,1)-x(3,end))-1,... % long vel
        abs(x(4:5,1)'-x(4:5,end)')-1,... % lat vel and yaw rate
        abs(x(6:7,1)'-x(6:7,end)')/0.330-1,...% front wheel speeds
        abs(x(8:9,1)'-x(8:9,end)')/0.330-1]';% rear wheel speeds

% c_if = abs(x(3,1)-x(3,end))-1;
% c_if = [];
c = [c; c_if ];

ceq = [ceq_dynamics' ceq_vehicle' ceq_initial-ceq_terminal];

end