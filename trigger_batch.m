% close;
clear;
clc;

%% Defines the folder name where the simulations results will be stored.
folder_name = 'test';
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

    optimalDesignFlag = false;
    
    problem = Initialize.fnInitMethod('trapezoid');
    problem = Initialize.fnInitVehicleSheet(problem,ii,setup_name4wm);
    problem = Initialize.fnInitTrack(problem, true, true);
    problem = Initialize.fnInitFunctionHandles(problem);
    problem = Initialize.fnInitBounds(problem,optimalDesignFlag);
    
    problem.options.optimalDesignFlag = optimalDesignFlag; 
    
    %% Get initial estimate
    % You need an initial guess for solving the BVP
    
    problem = PreProcessing.fnGetInitialEstimate(problem, 'PreSim');
    
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
    % options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
    %         'useparallel', true, ...
    %         'StepTolerance', 1e-5, ...
    %         'ConstraintTolerance', 1e-5, ...
    %         'OptimalityTolerance', 1e-5, ...
    %         'SpecifyConstraintGradient', false, ...
    %         'SpecifyObjectiveGradient', false, ...
    %         'MaxIterations', 1e4, ...
    %         'MaxFunctionEvaluations', 1e8, ...
    %         'FiniteDifferenceType','forward',...
    %         'FiniteDifferenceStepSize', 1e-5, ...
    %         'Algorithm','interior-point',...
    %         'ScaleProblem',true,...
    %         "EnableFeasibilityMode",true,...
    %         'PlotFcn', 'optimplotfvalconstr');
        
    problem.options.nlpOpt = options;    
    
    % Calls function that solves the NPL problem related to the MTM
    soln = Solver.MTM(problem);
    
    x_star = soln.grid.state;
    u_star = soln.grid.control;
    p_star = soln.grid.design;    
    Vehicle = problem.dsSystem.vd;
    Track = problem.dsSystem.td;
    
    %% getting additional outputs from optimal solution
    
    
    [~,~,~,O] = fnDynamicsVehicle([],x_star(3:end,:),u_star,Vehicle,p_star);
    
    
    KPIs = metrics_output(O,problem,soln);
    
    simResults(ii).solution = soln; 
    simResults(ii).outputs = O; 
    simResults(ii).metrics = KPIs; 

end

