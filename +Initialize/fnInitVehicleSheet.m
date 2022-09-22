function [problem] = fnInitVehicleSheet(problem,index,setup_name)

% this function loads the vehicle object to be used throughout the
% optimization

vehicle = Model.load_vehicle_parameters_from_sheet(index,setup_name);

problem.dsSystem.vd = vehicle;

end