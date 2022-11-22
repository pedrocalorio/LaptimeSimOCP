function [ceq,c] = fnCollectConstraints(defects,g,q,x,u,saving_constraints,saving_constraints_bound)

% collect deffects from collocation method to put in the ceq matrix of
% fmincon
ceq_dynamics = reshape(defects,numel(defects),1);
% ceq_dynamics = [];
% for i=0:size(defects,1)-1
%     ceq_dynamics = [ceq_dynamics, defects(i+1,:)];
% end

% collect the equality constraints on wheel load distribution
ceq_vehicle = reshape(q,numel(q),1);
% ceq_vehicle = [];
% for i=0:size(q,1)-1
%     ceq_vehicle = [ceq_vehicle, q(i+1,:)];
% end

% constraints so that the beggning and the end of the lap is the same
ceq_initial  = [x(1,1) x(2,1) x(3,1) x(4,1) x(5,1) x(8,1) x(9,1) x(10,1) x(11,1) u(1,1) u(2,1) u(3,1) u(4,1)];

ceq_terminal = [x(1,end) x(2,end) x(3,end) x(4,end) x(5,end) x(8,end) x(9,end) x(10,end) x(11,end) u(1,end) u(2,end) u(3,end) u(4,end)];

% inequality constraints on having always positve cornering stiffness and
% slip stiffness 
c = reshape(g,numel(g),1);
% c = [];

if saving_constraints_bound.fuel ~= -1
    c = [c; saving_constraints.fuel-saving_constraints_bound.fuel];
end

if saving_constraints_bound.tire_energy ~= -1
    c = [c; saving_constraints.tire_energy-saving_constraints_bound.tire_energy];
end


% ceq = [ceq_dynamics, ceq_vehicle, ceq_initial-ceq_terminal];

ceq = [ceq_dynamics' ceq_vehicle' ceq_initial-ceq_terminal];

end