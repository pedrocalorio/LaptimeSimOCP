% close;
clear;
clc;

%% Defines the folder name where the simulations results will be stored.

folder_name = 'Bahrein_v2';


% Folder path where the simulation results will be stored
results_folder = [pwd '\' 'SimResults'];

% Creates the folder to store results
mkdir(results_folder,folder_name);

% Gets the folder adress to later save plots and metrics table
base_folder_adress = [results_folder '\' folder_name];

%% Loads the speadsheet that contains vehicle parameters input
% Path of where the spreadsheet is. DO NOT CHANGE 

inputFolder = [pwd '\+Model'];

% Name of the spredsheet file. Change based on what is the vehicle parameters you want to load 
setup_name4wm = [inputFolder '\LMP2_OS.xlsx'];

% Gets the number of simulations based on the amount of different vehicles defined
number_of_sims = height(readtable(setup_name4wm,'Sheet','MassInertia'));

% Create a struct that will store simulation results and indexes from a batch run
simResults = struct;

%% Starts the simulation

% reads vehicle parameters from the spreadsheet
disp('Reading vehicle inputs...')
for ii = 1:number_of_sims

    %% Initialise the problem data structure
    
    % chooses what is going to be the integration scheme,
    % defines the number of discretization points
    %         if the simulation will optimize vehicle parameters
    %         if the simulation contains fuel saving constraints
    %         if the simulation contains tire energy saving constraints
    problem = Initialize.fnInitMethod('hermiteSimpson',600,false,false,false);

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
            'StepTolerance', 1e-6, ...
            'ConstraintTolerance', 1e-5, ...
            'OptimalityTolerance', 1e-6, ...
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

%         options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
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
%             'PlotFcn', 'optimplotfvalconstr');
        
    problem.options.nlpOpt = options;    
    
    % Calls function that solves the NPL problem related to the MTM
    disp("Starting the simulation number "+ii+'.')
    soln = Solver.MTM(problem);    


    simResults(ii).solution = soln; 

end

%% Post-processing loop
disp('Calculating post-processing toolbox...')
for ii = 1:number_of_sims
    % getting additional outputs from optimal solution

    soln = simResults(ii).solution;

    x_star  = soln.grid.state;
    u_star  = soln.grid.control;
    p_star  = soln.grid.design;    
    Vehicle = problem.dsSystem.vd;
    Track   = problem.dsSystem.td;
    
    
    [~,~,~,O,fuel_con] = Controller.fnDynamicsVehicle([],x_star(3:end,:),u_star,Vehicle,p_star,Track,x_star(1:2,:));
    
    
    KPIs = PostProcessing.metrics_output(O,problem,soln);
    
    simResults(ii).outputs = O; 
    simResults(ii).metrics = KPIs; 

end

%%

% saves the simulation results and all the variables in the workspace in a
% .mat file for loading it if necessary 

save(base_folder_adress+"\"+"sim.mat")