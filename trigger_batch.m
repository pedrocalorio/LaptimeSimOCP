% close;
clear;
clc;

%% Defines the folder name where the simulations results will be stored.
folder_name = 'fuji_lmp2_500pts_fuelsaving_v3_-25percent';
% folder_name = 'test_with_velocity';

% Folder path where the simulation results will be stored
results_folder = 'D:\dev\LaptimeSimOCP\SimResults';

% Creates the folder to store results
mkdir(results_folder,folder_name);

% Gets the folder adress to later save plots and metrics table
base_folder_adress = [results_folder '\' folder_name];

%% Loads the speadsheet that contains vehicle parameters input
% Path of where the spreadsheet is 
inputFolder = 'D:\dev\LaptimeSimOCP\+Model';

% Name of the spredsheet file
setup_name4wm = [inputFolder '\template.xlsx'];

% Gets the number of simulations based on the amount of different vehicles
% defined
number_of_sims = height(readtable(setup_name4wm,'Sheet','MassInertia'));

% Create a struct that will store simulation results and indexes from a
% batch run
simResults = struct;

%% Starts the simulation

% reads vehicle parameters from the spreadsheet
disp('Reading vehicle inputs...')
for ii = 1:number_of_sims

    %% Initialise the problem data structure
    
    problem = Initialize.fnInitMethod('trapezoid',false,true);
    problem = Initialize.fnInitVehicleSheet(problem,ii,setup_name4wm);
    problem = Initialize.fnInitTrack(problem, true, false);
    problem = Initialize.fnInitFunctionHandles(problem);
    problem = Initialize.fnInitBounds(problem);
    
    
    %% Get initial estimate
    % You need an initial guess for solving the BVP
    
    problem = PreProcessing.fnGetInitialEstimate(problem, 'Load');
    
    %% Set solver options and solve the Minimum Time Maneuvre 
    
    options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
            'useparallel', true, ...
            'StepTolerance', 1e-5, ...
            'ConstraintTolerance', 1e-5, ...
            'OptimalityTolerance', 1e-5, ...
            'SpecifyConstraintGradient', false, ...
            'SpecifyObjectiveGradient', false, ...
            'MaxIterations', 1e4, ...
            'MaxFunctionEvaluations', 1e8, ...
            'FiniteDifferenceType','forward',...
            'FiniteDifferenceStepSize', 1e-5, ...
            'Algorithm','interior-point',...
            'ScaleProblem',true,...
            "SubproblemAlgorithm","cg",...
            'PlotFcn', 'optimplotfvalconstr');
    % 
%     options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
%             'useparallel', true, ...
%             'StepTolerance', 1e-5, ...
%             'ConstraintTolerance', 1e-5, ...
%             'OptimalityTolerance', 1e-5, ...
%             'SpecifyConstraintGradient', false, ...
%             'SpecifyObjectiveGradient', false, ...
%             'MaxIterations', 1e4, ...
%             'MaxFunctionEvaluations', 1e8, ...
%             'FiniteDifferenceType','forward',...
%             'FiniteDifferenceStepSize', 1e-5, ...
%             'Algorithm','interior-point',...
%             'ScaleProblem',true,...
%             "EnableFeasibilityMode",true,...
%             'PlotFcn', 'optimplotfvalconstr');
        
    problem.options.nlpOpt = options;    
    
    % Calls function that solves the NPL problem related to the MTM
    soln = Solver.MTM(problem);
    


    simResults(ii).solution = soln; 

end
%% post processing loop
for ii = 1:number_of_sims
    % getting additional outputs from optimal solution

    soln = simResults(ii).solution;

    x_star  = soln.grid.state;
    u_star  = soln.grid.control;
    p_star  = soln.grid.design;    
    Vehicle = problem.dsSystem.vd;
    Track   = problem.dsSystem.td;
    
    
    [~,~,~,O,fuel_con] = fnDynamicsVehicle([],x_star(3:end,:),u_star,Vehicle,p_star,Track,x_star(1:2,:));
    
    
    KPIs = metrics_output(O,problem,soln);
    
    simResults(ii).outputs = O; 
    simResults(ii).metrics = KPIs; 

end

%%
save(base_folder_adress+"\"+folder_name)