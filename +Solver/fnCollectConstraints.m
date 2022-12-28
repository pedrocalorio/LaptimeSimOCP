function [ceq,c] = fnCollectConstraints(defects,g,x,u,...
    saving_constraints,saving_constraints_bound)

% collect deffects from collocation method to put in the ceq matrix of
% fmincon
ceq_dynamics = reshape(defects,numel(defects),1);


% constraints so that the beggning and the end of the lap is the same
ceq_initial  = [x(1,1) x(2,1) x(17,1) x(18,1) x(22,1) x(29,1) x(30,1) x(31,1) x(32,1) u(1,1) u(2,1) ];
ceq_terminal = [x(1,end) x(2,end) x(17,end) x(18,end) x(22,end) x(29,end) x(30,end) x(31,end) x(32,end) u(1,end) u(2,end) ];

% ceq_initial  = [x(:,1); u(:,1)]';
% ceq_terminal = [x(:,end); u(:,end)]';


% Inequality constraints on having always positve cornering stiffness and slip stiffness 
% c = reshape(g,numel(g),1);
% 
% if saving_constraints_bound.fuel ~= -1
%     c = [c; saving_constraints.fuel-saving_constraints_bound.fuel];
% end
% 
% if saving_constraints_bound.tire_energy ~= -1
%     c = [c; (saving_constraints.tire_energy-saving_constraints_bound.tire_energy)];
% end
c=[];



ceq = [ceq_dynamics' ceq_initial-ceq_terminal];

end