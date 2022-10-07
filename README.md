# LaptimeSim with OCP

This is the main repository for the laptime simulation using optimal control created in MATLAB. In this file, I'm going to explain how you can run a simulation successfully 
and what are the features you can control

## How to run a simulation?

You should be using the script "trigger_batch" to call the functions that is going 
to be solving the problem

### Loading vehicle parameters 

The vehicle parameter inputs are loaded from an excel sheet. You can find one sheet 
template in "+Model/template.xlsx". In case you want to do a batch run with multiple
vehicles using different settings, you just add another entry in the sheet in a new
row.

When creating a new excel file, please put them inside the folder "+Model" because that
is the folder the software will look for to load the vehicle parameter inputs.

Two important comments:

1. Tire model is not loaded inside the sheet, instead it is always loaded by referring to 
a script which loads all the tire model coefficients. You can check how that is working
in the LMPTire_Front_26psi script.
2. The engine map is loaded from an excel spreadsheet to and it is also inside the 
"+Model" folder.

### Loading the track

The track is generated when calling the function fnInitTrack. For 
this function you pass 3 arguments: the first is the problem struct which is going 
to store all the information needed to solve the simulation. The second one is a 
boolean flag if you want to plot the track layout, in case you want to see the plots for the
track this should be set to "true". 

Lastly, and more importantly, we have another boolean getting the
information whether the track is going to be re-created from telemetry data (set it to false) or
if you want to construct an oval track using mathematics and geometry (set to true). 

In case the option of re-constructing the track from data is chosen, it is needed 
to go inside the function fnInitTrack, go on line 81 that will load the data and set
the name of the .csv file that contains this data. It is important that the format
of the sheet follows this order of channels:

Time[s] || Distance[m] || Speed[km/h] || Lateral Acceleration [g]

The function generate_track_from_telemetry is responsible to do all the processing
in which the track curvature, x, y and yaw_angle variables are going to be calculated
and stored in the Track struct. Those variables are very important for the solver.

### Initial Guess

For this kind of simulation, an initial guess is required and this is defined in the function
fnGetInitialEstimate. 

It can be defined as 'PreSim' in which the software is going to calculate
a relatively poor initial guess and the optimizer might need more time until convergence. The
other option is to use 'Load', in which the optimizer will use a previous simulation
as initial guess. This can really speed up the process and could be useful if you already have 
solved the laptime for that specific track and just want to simulate it again for a different
vehicle setup.

In general, if you are simulating a track for the first time, it is recommended to use the PreSim
Otherwise, you should use Load options because it will significantly reduce time.

### Define what type of laptime simulation to run

This program solves the minimum time of a given vehicle on a given track. However,
there are different ways under different conditions that you can achive that. So,
there are 3 flags that you need to define which are going set the way in which the 
problem is going to be solved, and they are:
   1. If you would like to optimize vehicle parameters, and for the moment those are
the parameters we have implemented:
      a. Weight distribution
      b. Downforce distribution
      c. Lateral load transfer distribution
      d. Brake torque distribution
   2. If the simulation has fuel saving constraints
   3. If the simulation has tire energy saving constraints

In case one of them is set to true, the admissible variables are defined inside the
function fnInitBounds

## Considering vehicle stability in the problem

When running simulations with vehicle parameter optimization, you can include
constraints on the stability of the vehicle, so that the optimal solution doesn't
yield in an undrivable vehicle. Those constraints are included by
simply calculating the eigenvalues of an estimated simple bicycle model at each
iteration and making sure that there are no right-side poles. This idea is better
explained in my thesis, in case you are interested.

That being said, in case you want to include them, it is just a matter of uncommenting line 290 
of the fnDynamicsVehicle function. I'm still thinking of better calculations for 
considering stability in this problem, for the moment is very hard-coded but this 
is something that I would like to change.

Once you have the vehicle parameters ready in the spreadsheet, the track to be loaded 
is defined, the kind of laptime simulation is also set, then you are ready to finally
run the simulation and those should be the steps:

1. Define the name of the folder in "folder_name" that is going to be created
   with the simulation results
2. Select the name of the sheet that is going to load the vehicle parameters input 
   for the simulation and put it in the "setup_name4wm" variable 

After the simulation ends, you can call the ploting_script to display some of the results.

But in general, the main output of the simulation is stored in the simResults struct.
In there you have a field for metrics, which are some pre-calculated metrics around 
this optimal lap like steering integral, laptime and tire energies. You also have 
a matrix with some pre-defined outputs like normal loads, lateral force, yaw moments, etc.
Finally you also have a solution field that is going to store every information on 
the status of the optimization and the settings of the simulation that you have set
previously. 

In the case of a batch-run, there will be a list of those parameters for each run. 
For example, if 4 simulations are defined to run in the sheet, the simResults struct 
will have those 3 fields, but one entry for each setup, so if you enter in the metrics 
field, there will be array with 4 numbers, each one corresponding to each run.
