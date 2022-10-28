function [problem] = fnPreSimulate(problem)

td = problem.dsSystem.td;
vd = problem.dsSystem.vd;

curv = interp1(td.distance,td.curv,td.sLap,'spline');

% Vx0   = min(0.08.*td.sLap, (0.5*9.81./abs(td.curv)).^0.5);
Vx0 = min(180/3.6, (1.5*9.81./abs(curv)).^0.5);
% Vx0   = zeros(1,length(td.curv));
beta0 = zeros(1,length(Vx0));
r0    = Vx0.*curv;
omega10 = Vx0./vd.tire_1.radius;
omega20 = Vx0./vd.tire_2.radius;
omega30 = Vx0./vd.tire_3.radius;
omega40 = Vx0./vd.tire_4.radius;

n0    = zeros(1,length(Vx0));
zeta0 = zeros(1,length(Vx0));


delta0 = 0.*curv;
deltaRate0 = gradient(delta0);
% delta0 = zeros(1,length(Vx0));
tau   = 0.0*ones(1,length(Vx0));
tauRate   = gradient(tau);
LongLT0   = zeros(1,length(Vx0)); 
LatLT0    = zeros(1,length(Vx0)); 

problem.guess.design = [problem.dsSystem.vd.chassis.frontMassDist;
    problem.dsSystem.vd.aero.aeroBalance;
    problem.dsSystem.vd.susp.LatLTD;
    problem.dsSystem.vd.brakes.bias];

problem.guess.state   	= [problem.dsSystem.PreSim.scale.*n0;
                           problem.dsSystem.PreSim.scale.*zeta0; 
                           problem.dsSystem.PreSim.scale.*Vx0;
                           problem.dsSystem.PreSim.scale.*beta0;
                           problem.dsSystem.PreSim.scale.*r0;
                           problem.dsSystem.PreSim.scale.*omega10;
                           problem.dsSystem.PreSim.scale.*omega20;
                           problem.dsSystem.PreSim.scale.*omega30;
                           problem.dsSystem.PreSim.scale.*omega40;
                           problem.dsSystem.PreSim.scale.*delta0;
                           problem.dsSystem.PreSim.scale.*tau];
                       
problem.guess.control   = [problem.dsSystem.PreSim.scale.*deltaRate0;
                           problem.dsSystem.PreSim.scale.*tauRate;
                           problem.dsSystem.PreSim.scale.*LongLT0
                           problem.dsSystem.PreSim.scale.*LatLT0];

end