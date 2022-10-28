function vehicle = load_vehicle_parameters()

vehicle = struct();

% Chassis

vehicle.chassis.Ms              = 1030;
vehicle.chassis.frontMassDist   = 47.0/100;
vehicle.chassis.wheelbase       = 3.05;
% vehicle.chassis.a = vehicle.chassis.wheelbase * (1-vehicle.chassis.frontMassDist);
% vehicle.chassis.b = vehicle.chassis.wheelbase - vehicle.chassis.a;
vehicle.chassis.frontTrack      = 1.560;
vehicle.chassis.rearTrack       = 1.550;
vehicle.chassis.Iz              = 1250;
vehicle.chassis.CoG             = 0.350;

% Aero 

vehicle.aero.liftCoeff          = 4.0;
vehicle.aero.dragCoeff          = 1.000;
vehicle.aero.frontalArea        = 1.0;
vehicle.aero.airDensity         = 1.226;
vehicle.aero.aeroBalance        = 44/100;

% Susp

vehicle.susp.LatLTD  = 44/100; 

% Steering

vehicle.steeringratio = 13.71;

% % Non-suspeded mass properties
% 
% vehicle.mu1 = 21;
% vehicle.mu2 = 21;
% vehicle.mu3 = 24;
% vehicle.mu4 = 24;

vehicle.tire_1.inertia = 3;
vehicle.tire_2.inertia = 3;
vehicle.tire_3.inertia = 3;
vehicle.tire_4.inertia = 3;

% % Tire properties
% 
vehicle.tire_1.radius = 0.330;
vehicle.tire_2.radius = 0.330;
vehicle.tire_3.radius = 0.330;
vehicle.tire_4.radius = 0.330;

% vehicle.tire_1.MF = mfeval.readTIR('637319_ISO_FRONT.TIR');
% vehicle.tire_2.MF = mfeval.readTIR('637319_ISO_FRONT.TIR');
% vehicle.tire_3.MF = mfeval.readTIR('637326_ISO_REAR.TIR');
% vehicle.tire_4.MF = mfeval.readTIR('637326_ISO_REAR.TIR');

vehicle.tire_1.MF = Model.LMPTire_Front_26psi_v2();
vehicle.tire_2.MF = Model.LMPTire_Front_26psi_v2();
vehicle.tire_3.MF = Model.LMPTire_Rear_26psi_v2();
vehicle.tire_4.MF = Model.LMPTire_Rear_26psi_v2();

% Diff

vehicle.differential.To = 1*50;
vehicle.differential.G = 1.8;

% Engine properties

EnginePowerMap = readtable("D:\dev\LaptimeSimOCP\+Model\LMP2 Engine Power Map.xlsx");

% going from hp to watts the multiplicy constant is 745.7
vehicle.engine.maximum_power = 550*745.7;
vehicle.engine.power_drag = 0*745.7;

vehicle.engine.RPM = EnginePowerMap.RPM;
vehicle.engine.Power_R = EnginePowerMap.Power_R;

vehicle.engine.Throttle = EnginePowerMap.Throttle;
vehicle.engine.Power_T = EnginePowerMap.Power_T;

% vehicle.engine.Gr = 1.05*2.9; % gear ratio which is constant for the moment
vehicle.engine.Gr = 1.30*2.9; % gear ratio which is constant for the moment

% Brake properties

vehicle.brakes.bias = 0.70;
vehicle.brakes.maximum_torque = 5000;


end

