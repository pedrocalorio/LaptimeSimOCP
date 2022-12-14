function vehicle = load_vehicle_parameters_from_sheet(index,setup_name)

vehicle = struct();

[MassInertia, Dimensions, Kinematics, Aerodynamics, Suspension,...
    TireFront, TireRear, Brake, Differential, Engine] = Model.load_vehicle_table(setup_name);

% Chassis
vehicle.chassis.Ms              = MassInertia.Ms(index);
vehicle.chassis.frontMassDist   = MassInertia.weight_d(index);
vehicle.chassis.wheelbase       = Dimensions.wheelbase(index);
vehicle.chassis.frontTrack      = Dimensions.front_track(index);
vehicle.chassis.rearTrack       = Dimensions.rear_track(index);
vehicle.chassis.Iz              = MassInertia.Iz(index);
vehicle.chassis.CoG             = MassInertia.h_cg(index);

% Aero 
vehicle.aero.liftCoeff          = Aerodynamics.C_down(index);
vehicle.aero.dragCoeff          = Aerodynamics.C_drag(index);
vehicle.aero.frontalArea        = Aerodynamics.frontal_area(index);
vehicle.aero.airDensity         = Aerodynamics.air_dens(index);
vehicle.aero.aeroBalance        = Aerodynamics.aerobalance(index);

% Susp
vehicle.susp.LatLTD  = Suspension.LLTD(index); 

% Steering
vehicle.steeringratio = Kinematics.i_s(index);

% Camber
vehicle.camber_s_1 = Kinematics.static_camber_1(index);
vehicle.camber_s_2 = Kinematics.static_camber_2(index);
vehicle.camber_s_3 = Kinematics.static_camber_3(index);
vehicle.camber_s_4 = Kinematics.static_camber_4(index);

% here I say that the wheel rotational inertia is a tire property, but it
% is not true
vehicle.tire_1.inertia = MassInertia.I_yp_1(index);
vehicle.tire_2.inertia = MassInertia.I_yp_2(index);
vehicle.tire_3.inertia = MassInertia.I_yp_3(index);
vehicle.tire_4.inertia = MassInertia.I_yp_4(index);

% % Tire properties
% 
vehicle.tire_1.radius = TireFront.radius(index);
vehicle.tire_2.radius = TireFront.radius(index);
vehicle.tire_3.radius = TireRear.radius(index);
vehicle.tire_4.radius = TireRear.radius(index);

vehicle.tire_1.rr = TireFront.fr(index);
vehicle.tire_2.rr = TireFront.fr(index);
vehicle.tire_3.rr = TireRear.fr(index);
vehicle.tire_4.rr = TireRear.fr(index);

% vehicle.tire_1.MF = mfeval.readTIR('637319_ISO_FRONT.TIR');
% vehicle.tire_2.MF = mfeval.readTIR('637319_ISO_FRONT.TIR');
% vehicle.tire_3.MF = mfeval.readTIR('637326_ISO_REAR.TIR');
% vehicle.tire_4.MF = mfeval.readTIR('637326_ISO_REAR.TIR');

vehicle.tire_1.MF = Model.LMPTire_Front_26psi_v2();
vehicle.tire_2.MF = Model.LMPTire_Front_26psi_v2();
vehicle.tire_3.MF = Model.LMPTire_Rear_26psi_v2();
vehicle.tire_4.MF = Model.LMPTire_Rear_26psi_v2();

% Diff

vehicle.differential.To = Differential.preload(index);
vehicle.differential.G = Differential.sens(index);

% Engine properties

EnginePowerMap = readtable("+Model\"+string(Engine.map(index)));

% going from hp to watts the multiplicy constant is 745.7
vehicle.engine.maximum_power = Engine.maximum_power(index)*745.7;
vehicle.engine.power_drag = Engine.power_drag(index)*745.7;

vehicle.engine.RPM = EnginePowerMap.RPM;
vehicle.engine.Power_R = EnginePowerMap.Power_R;

vehicle.engine.Throttle = EnginePowerMap.Throttle;
vehicle.engine.Power_T = EnginePowerMap.Power_T;

vehicle.engine.Gr = readmatrix("+Model\"+string(Engine.gear_ratio(index))); % gear ratio which is constant for the moment

% Brake properties

vehicle.brakes.bias = Brake.bias(index);
vehicle.brakes.maximum_torque = Brake.maximum_torque(index);


end

