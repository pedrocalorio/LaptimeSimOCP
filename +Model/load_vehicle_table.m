function [MassInertia, Dimensions, Kinematics, Aero, Susp,...
    TireFront, TireRear, Brake, Diff, Engine] = load_vehicle_table(setup_name)

%INPUT----
%- setup_name ( string, [-])            String name that contains the name
%                                       of the speadsheet in which the vehicle 
%                                       parameters are stored.


%OUTPUT----
%- MassInertia ( table, [-])                Loads sprung mass; sprung mass
%  inertia in x,y and z; unsprung mass of 4 corners; wheel rotational inertia
%  of the 4 wheels; CG height of the sprung mass; weight distribution and
%  constant of gravity
%
%- Dimensions ( table, [-])                 Loads wheelbase; wheel tracks and
%  unloaded radius 
%
%- Kinematics ( table, [-])                 Loads the front and side view
%  instant center angles (for jacking forces); camber and toe angles along
%  with the steering ratio.
%
%- Aero ( table, [-])                       Loads aerodynamics parameters such
%  as drag and downforce coefficient; x,y and z cp coordinate 
%
%- Tire [front and rear] ( table, [-])      Loads tire properties like
%  lateral and longitudinal relaxation length; rolling resistance coefficients;
%  tire model string; (vertical) damping and stiffness; pressure;
%  freeLength and solidLength
%
%- Damper [corners and heave] ( table, [-])   Loads the damping
%  properties such as the stiffness; freeLength; solidLength; static pressure;
%  shaft diameter and motion ratio 
%
%- Bumpstop, Spring and ARB [corners and heave] ( table, [-])   Loads the
%  elastic element properties (the same for bumpstop and spring) such as
%  stiffness; freeLength; solidLength; gaps and motion ratio 

% Load mass and inertia properties from the table
MassInertia = readtable(setup_name,'Sheet','MassInertia');

% Load dimensions properties from the table
Dimensions = readtable(setup_name,'Sheet','Dimensions');

% Load Kinematics properties from the table
Kinematics = readtable(setup_name,'Sheet','Kinematics');

% Load Aero properties from the table
Aero = readtable(setup_name,'Sheet','Aero');

% Load Susp properties from the table
Susp = readtable(setup_name,'Sheet','Susp');

% Loads tire related properties from the table
TireFront = readtable(setup_name,'Sheet','TireFront');
TireRear = readtable(setup_name,'Sheet','TireRear');

% Loads tire related properties from the table
Brake = readtable(setup_name,'Sheet','Brake');

Diff = readtable(setup_name,'Sheet','Diff');

Engine = readtable(setup_name,'Sheet','Engine');