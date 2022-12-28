function [problem] = fnPreSimulate(problem)

td = problem.dsSystem.td;
vd = problem.dsSystem.vd;

% curv = interp1(td.distance,td.curv,td.sLap,'spline');
curv = interp1(td.sLap,td.curv,td.sLap,'spline');

% Vx0   = min(0.08.*td.sLap, (0.5*9.81./abs(td.curv)).^0.5);
Vx0 = min(150/3.6, (1.0*9.81./abs(curv)).^0.5);
% Vx0   = zeros(1,length(td.curv));
beta0 = zeros(1,length(Vx0));
r0    = Vx0.*curv;
omega10 = Vx0./vd.tire_1.radius;
omega20 = Vx0./vd.tire_2.radius;
omega30 = Vx0./vd.tire_3.radius;
omega40 = Vx0./vd.tire_4.radius;

psi0 = cumtrapz(td.sLap,r0)/td.sLap(end);
long_pos = cumtrapz(td.sLap,Vx0)/td.sLap(end);
omega1_integral = cumtrapz(td.sLap,omega10)/td.sLap(end);
omega2_integral = cumtrapz(td.sLap,omega20)/td.sLap(end);
omega3_integral = cumtrapz(td.sLap,omega30)/td.sLap(end);
omega4_integral = cumtrapz(td.sLap,omega40)/td.sLap(end);



w0 = zeros(1,length(Vx0));
wx0 = zeros(1,length(Vx0));
wy0 = zeros(1,length(Vx0));

z1_dot0 = zeros(1,length(Vx0));
z2_dot0 = zeros(1,length(Vx0));
z3_dot0 = zeros(1,length(Vx0));
z4_dot0 = zeros(1,length(Vx0));

n0    = zeros(1,length(Vx0));
zeta0 = zeros(1,length(Vx0));

delta0 = 0.1*curv;
deltaRate0 = gradient(delta0);
% delta0 = zeros(1,length(Vx0));
tau   = 0.5*ones(1,length(Vx0));
tauRate   = gradient(tau);

% xtirf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
% xtilf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
% xtilr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
% xtirr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
% 
% xsirf=(m*g*b)/(2*(a+b)*ksrf);
% xsilf=(m*g*b)/(2*(a+b)*kslf);
% xsilr=(m*g*a)/(2*(a+b)*kslr);
% xsirr=(m*g*a)/(2*(a+b)*ksrr);

xtirf = ones(1,length(Vx0))*0.01; 
xtilf = xtirf;
xtilr = ones(1,length(Vx0))*0.015; 
xtirr = xtilr;

xsirf = ones(1,length(Vx0))*0.02; 
xsilf = xsirf;
xsilr = ones(1,length(Vx0))*0.025; 
xsirr = xsilr;

problem.guess.design = [problem.dsSystem.vd.chassis.frontMassDist;
    problem.dsSystem.vd.aero.aeroBalance;
    problem.dsSystem.vd.susp.LatLTD;
    problem.dsSystem.vd.brakes.bias];

problem.guess.state   	= [problem.dsSystem.PreSim.scale.*n0;
                           problem.dsSystem.PreSim.scale.*zeta0; 
                           problem.dsSystem.PreSim.scale.*long_pos;
                           zeros(2,length(Vx0));
                           problem.dsSystem.PreSim.scale.*0.05*curv;
                           zeros(1,length(Vx0));
                           problem.dsSystem.PreSim.scale.*psi0;
                           problem.dsSystem.PreSim.scale.*-xsilf;
                           problem.dsSystem.PreSim.scale.*-xsirf;
                           problem.dsSystem.PreSim.scale.*-xsilr;
                           problem.dsSystem.PreSim.scale.*-xsirr;
                           omega1_integral;
                           omega2_integral;
                           omega3_integral;
                           omega4_integral;
                           problem.dsSystem.PreSim.scale.*Vx0;
                           problem.dsSystem.PreSim.scale.*beta0;
                           problem.dsSystem.PreSim.scale.*w0;
                           problem.dsSystem.PreSim.scale.*wx0;
                           problem.dsSystem.PreSim.scale.*wy0;
                           problem.dsSystem.PreSim.scale.*r0;
                           problem.dsSystem.PreSim.scale.*z1_dot0;
                           problem.dsSystem.PreSim.scale.*z2_dot0;
                           problem.dsSystem.PreSim.scale.*z3_dot0;
                           problem.dsSystem.PreSim.scale.*z4_dot0;
                           problem.dsSystem.PreSim.scale.*omega10;
                           problem.dsSystem.PreSim.scale.*omega20;
                           problem.dsSystem.PreSim.scale.*omega30;
                           problem.dsSystem.PreSim.scale.*omega40;
                           problem.dsSystem.PreSim.scale.*delta0;
                           problem.dsSystem.PreSim.scale.*tau;
                           problem.dsSystem.PreSim.scale.*xsilf;
                           problem.dsSystem.PreSim.scale.*xsirf;
                           problem.dsSystem.PreSim.scale.*xsilr;
                           problem.dsSystem.PreSim.scale.*xsirr;
                           problem.dsSystem.PreSim.scale.*xtilf;
                           problem.dsSystem.PreSim.scale.*xtirf;
                           problem.dsSystem.PreSim.scale.*xtilr;
                           problem.dsSystem.PreSim.scale.*xtirr];
                       
problem.guess.control   = [problem.dsSystem.PreSim.scale.*deltaRate0;
                           problem.dsSystem.PreSim.scale.*tauRate];

end