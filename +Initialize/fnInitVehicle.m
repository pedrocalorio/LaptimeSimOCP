function [problem] = fnInitVehicle(problem)

% this function loads the vehicle object to be used throughout the
% optimization

vehicle = Model.load_vehicle_parameters();

problem.dsSystem.vd = vehicle;

end