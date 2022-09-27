function [problem] = fnGetInitialEstimate(problem, sMethod)

switch sMethod
    case 'Load'
        %
%         customGuess = load([cd, '\+SimResults\girish tire model\BAR0_hs_300pts_no_engine_model_OD.mat']);
%         customGuess = load([cd, '\+SimResults\BAR_optimization_trap_1050pts_no_engine_model_NO_opt_des_base_v2_up.mat']);
        customGuess = load([cd, '\SimResults\fuji_lmp2_500pts_fuel_saving_v2\fuji_lmp2_500pts_fuel_saving_v2.mat']);

        sLap = customGuess.soln.grid.sLap;

        sLap_new = problem.dsSystem.td.sLap;

        scale = 1;

        if problem.options.optimalDesignFlag == 1
            problem.guess.design = [customGuess.soln.grid.design(1);
                    customGuess.soln.grid.design(2);
                    customGuess.soln.grid.design(3);
                    customGuess.soln.grid.design(4)];
        else
            problem.guess.design = [];
        end

       
         
        problem.guess.state   	= scale.* [interp1(sLap,customGuess.soln.grid.state(1,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(2,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(3,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(4,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(5,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(6,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(7,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(8,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(9,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(10,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.state(11,:),sLap_new)];
                               
        problem.guess.control   =scale.* [interp1(sLap,customGuess.soln.grid.control(1,:),sLap_new);                                   
                                   interp1(sLap,customGuess.soln.grid.control(2,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.control(3,:),sLap_new);
                                   interp1(sLap,customGuess.soln.grid.control(4,:),sLap_new)];

    case 'LoadData'
        customGuess = load([cd, '/#65_QU_RUN001-WS_21_04_17_BARCELONE.mat']);

        scale = 1;

        sLap = linspace(0,customGuess.Alias_Lap_Distance.Value(end)-customGuess.Alias_Lap_Distance.Value(1),length(customGuess.Alias_Throttle_Pos.Value));

        sLap_new = problem.dsSystem.td.sLap;

        factor = sLap(end)/sLap_new(end);
        sLap = sLap/factor;

        problem.guess.design = [problem.dsSystem.vd.chassis.frontMassDist;
            problem.dsSystem.vd.aero.aeroBalance;
            problem.dsSystem.vd.susp.LongLTD;
            problem.dsSystem.vd.susp.LatLTD;
            problem.dsSystem.vd.brakes.bias];

        tau = customGuess.Alias_Throttle_Pos.Value/1e2;
        tauneg = -customGuess.Math_Brake_Pres_Total.Value / max(customGuess.Math_Brake_Pres_Total.Value);

        for i=1:length(tau)
            if tau(i)<5/100
                tau(i) = tauneg(i);
            else
                tau(i) = tau(i);
            end
        end

        problem.guess.state   	= scale.* [zeros(1,length(sLap_new));
                           zeros(1,length(sLap_new));
                           interp1(sLap,customGuess.Corr_Speed.Value/3.6,sLap_new);
                           zeros(1,length(sLap_new));
                           interp1(sLap,customGuess.Alias_Gyro_Yaw_Velocity.Value/3.6,sLap_new);
                           interp1(sLap,customGuess.Corr_Speed.Value/3.6,sLap_new)./problem.dsSystem.vd.tire_1.radius;
                           interp1(sLap,customGuess.Corr_Speed.Value/3.6,sLap_new)./problem.dsSystem.vd.tire_2.radius;
                           interp1(sLap,customGuess.Corr_Speed.Value/3.6,sLap_new)./problem.dsSystem.vd.tire_3.radius;
                           interp1(sLap,customGuess.Corr_Speed.Value/3.6,sLap_new)./problem.dsSystem.vd.tire_4.radius;
                           interp1(sLap,customGuess.Alias_Steering_Wheel_Angle.Value/57.3,sLap_new);
                           interp1(sLap,tau,sLap_new)];

        problem.guess.control   =scale.* [gradient(interp1(sLap,customGuess.Alias_Steering_Wheel_Angle.Value/57.3,sLap_new));                                   
                           gradient(interp1(sLap,tau,sLap_new));
                           zeros(1,length(sLap_new));
                           zeros(1,length(sLap_new))];
            
            
        
    case 'PreSim'
        problem.dsSystem.PreSim.enabled = true;
        problem.dsSystem.PreSim.scale = 1;
        problem = PreProcessing.fnPreSimulate(problem);
    otherwise
        error('Select a method for finding the initial estimate');
end

end