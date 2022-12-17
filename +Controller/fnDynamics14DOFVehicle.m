function [dx,g_ineq,saving_constraints] = fnDynamics14DOFVehicle(x,u_c,p,Vehicle)

addpath('C:/dev/libraries/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

% Distance Step
ds = (Track.sLap(end) - Track.sLap(1)) / (length(x(1,:)) - 1);

% Gravity constant
g = 9.81;

%% Vehicle states
% heave = x(1); % longitudinal displacement
% y = x(2); % lateral displacement
heave = x(3,:); % heave displacement
phi   = x(4,:); %roll angle
theta = x(5,:); %pitch angle
psi   = x(6,:); %yaw angle
z1    = x(7,:); % fl unsprung mass displacement
z2    = x(8,:); % fr unsprung mass displacement
z3    = x(9,:); % rl unsprung mass displacement
z4    = x(10,:); % rr unsprung mass displacement

%%%%%%%%%%%%%%%%%%%%%%%% I don't need the wheel angular position

u  = x(15,:);
v  = x(16,:);
w  = x(17,:);
wx = x(18,:);
wy = x(19,:);
wz = x(20,:);
z1_dot = x(21,:); % fl unsprung mass displacement rate
z2_dot = x(22,:); % fr unsprung mass displacement rate
z3_dot = x(23,:); % rl unsprung mass displacement rate
z4_dot = x(24,:); % rr unsprung mass displacement rate
omega1 = x(25,:); 
omega2 = x(26,:); 
omega3 = x(27,:); 
omega4 = x(28,:); 
steer = x(29,:);
tau   = x(30,:);

%% Control inputs
steerRate = u_c(1,:);
tauRate = u_c(2,:);

%% Calculating the cumulative laptime

Sf = (1 - n.*kappa)./(vx.*cos(xi)-vy.*sin(xi));
cumLaptime = ds*Sf*ones(length(n),1);

%% Vehicle parameters
% Aero
liftCoeff   = Vehicle.aero.liftCoeff;
dragCoeff   = Vehicle.aero.dragCoeff;
frontalArea = Vehicle.aero.frontalArea;
airDens     = Vehicle.aero.airDensity;
aeroBalance = p(2);

% Susp
kslf=   160e3; %front left suspension stiffness (N/m)
ksrf=   160e3;
kslr=   95e3;  %rear left suspension stiffness (N/m)
ksrr=   95e3;
% bslf=   7000;  %front left suspension damping coefficient (Ns/m)
% bsrf=   7000;
bslf=   p(3);  %front left suspension damping coefficient (Ns/m)
bsrf=   p(3);
bslr=   7000;  %rear left suspension damping coefficient (Ns/m)
bsrr=   7000;  
hrcf=   0.2;   %front roll center distance below sprung mass c.g. (m)
hrcr=   0.15;  %rear roll center distance below sprung mass c.g. (m)

% Steering
steeringRatio = Vehicle.steeringratio;

% Chassis
m=Vehicle.chassis.Ms; % Sprung mass (kg)
Jx=900; % Sprung mass roll inertia (kg.m^2)
Jy=2000; % Sprung mass pitch inertia (kg.m^2)
Jz=Vehicle.chassis.Iz; % Sprung mass yaw inertia (kg.m^2)
wheelbase   = Vehicle.chassis.wheelbase;
weight_dist = p(1);
a           = wheelbase * (1-weight_dist);
b           = wheelbase * weight_dist;
h           = Vehicle.chassis.CoG;  % Sprung mass c.g. height (m)
cf = Vehicle.chassis.frontTrack;
cr = Vehicle.chassis.rearTrack; % front/rear track width (m)
muf=55;    %front unsprung mass (kg)
mur=55;    %rear unsprung mass (kg)

% Tires
r01 = Vehicle.tire_1.radius;
r02 = Vehicle.tire_2.radius;
r03 = Vehicle.tire_3.radius;
r04 = Vehicle.tire_4.radius;
Iwheel_1 = Vehicle.tire_1.inertia;
Iwheel_2 = Vehicle.tire_2.inertia;  
Iwheel_3 = Vehicle.tire_3.inertia;
Iwheel_4 = Vehicle.tire_4.inertia; 
ktf=300e3;  %front tire stiffness (N/m)
ktr=300e3;  %rear tire stiffness (N/m)
MF_f = Vehicle.tire_1.MF;
MF_r = Vehicle.tire_3.MF;   

%% Calculations

% steering angle at the tire
delta = steer/steeringRatio;

% the initial tire compression xtif
xtirf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
xtilf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
xtilr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
xtirr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
xtlf=xtilf-z1;
xtrf=xtirf-z2;
xtlr=xtilr-z3;
xtrr=xtirr-z4;

%the initial spring compression xsif 
xsirf=(m*g*b)/(2*(a+b)*ksrf);
xsilf=(m*g*b)/(2*(a+b)*kslf);
xsilr=(m*g*a)/(2*(a+b)*kslr);
xsirr=(m*g*a)/(2*(a+b)*ksrr);
xsrf=xsirf + heave + a*sin(theta) + cf/2*sin(phi);
xslf=xsilf + heave + a*sin(theta) - cf/2*sin(phi);
xslr=xsilr + heave - b*sin(theta) + cr/2*sin(phi);
xsrr=xsirr + heave - b*sin(theta) - cr/2*sin(phi);

% the initial length of the strut
% lsirf=h-(r0-xtirf);
% lsilf=h-(r0-xtilf);
% lsilr=h-(r0-xtilr);
% lsirr=h-(r0-xtirr);
lsrf=h-(r01-xtrf)+heave;
lslf=h-(r02-xtlf)+heave;
lslr=h-(r03-xtlr)+heave;
lsrr=h-(r04-xtrr)+heave;

% the transforming matrix 
Mrf=[0 0 cf/2;0 0 a;-cf/2 -a 0];
Mlf=[0 0 -cf/2;0 0 a;cf/2 -a 0];
Mlr=[0 0 -cr/2;0 0 -b;cr/2 b 0];
Mrr=[0 0 cr/2;0 0 -b;-cr/2 b 0];
% initial velocity of strut mounting point in x y z by transforming the c.g. velocities
wm=[wx;wy;wz];
vm=[u;v;w];

% vou precisar fazer um for loop
for i=1:length(x(1,:))
    Vsrf(:,1)=Mrf*wm(:,i)+vm(:,i); % front right mounting point x y z velocity in coordinate frame 1
    Vslf(:,1)=Mlf*wm(:,i)+vm(:,i); % front left
    Vslr(:,1)=Mlr*wm(:,i)+vm(:,i); % rear left
    Vsrr(:,1)=Mrr*wm(:,i)+vm(:,i); % rear right
end
% Vsrf=Mrf*wm+vm; % front right mounting point x y z velocity in coordinate frame 1
% Vslf=Mlf*wm+vm; % front left
% Vslr=Mlr*wm+vm; % rear left
% Vsrr=Mrr*wm+vm; % rear right

% initial unsprung mass velocity
uurf=Vsrf(1,:)-lsrf.*wy; 
uulf=Vslf(1,:)-lslf.*wy; 
uulr=Vslr(1,:)-lslr.*wy; 
uurr=Vsrr(1,:)-lsrr.*wy; 
vurf=Vsrf(2,:)+lsrf.*wx; 
vulf=Vslf(2,:)+lslf.*wx; 
vulr=Vslr(2,:)+lslr.*wx; 
vurr=Vsrr(2,:)+lsrr.*wx;
wulf=w;
wurf=w;
wulr=w;
wurr=w;
% % initial strut acceration
dxslf=-Vslf(3,:)+wulf;
dxsrf=-Vsrf(3,:)+wurf;
dxslr=-Vslr(3,:)+wulr;
dxsrr=-Vsrr(3,:)+wurr;
% initial unsprung mass acceration


% wheel rotational velocity
% wlf=u/(r0-xtlf);
% wrf=u/(r0-xtrf);
% wlr=u/(r0-xtlr);
% wrr=u/(r0-xtrr);


% the instantaneous tire radius
% % to account for the wheel lift-off, when the tire radial compression becomes less than zero, Rij=r0;

Rrf=(r01-xtrf)./(cos(theta).*cos(phi));

Rlf=(r02-xtlf)./(cos(theta).*cos(phi));

Rlr=(r03-xtlr)./(cos(theta).*cos(phi));

Rrr=(r04-xtrr)./(cos(theta).*cos(phi));

% the longitudinal and lateral velocities at the tire contact patch in coordinate frame 2
ugrf=cos(theta).*(uurf-wy.*Rrf)+sin(theta).*(wurf.*cos(phi)+sin(phi).*(wx.*Rrf+vurf));
uglf=cos(theta).*(uulf-wy.*Rlf)+sin(theta).*(wulf.*cos(phi)+sin(phi).*(wx.*Rlf+vulf));
uglr=cos(theta).*(uulr-wy.*Rlr)+sin(theta).*(wulr.*cos(phi)+sin(phi).*(wx.*Rlr+vulr));
ugrr=cos(theta).*(uurr-wy.*Rrr)+sin(theta).*(wurr.*cos(phi)+sin(phi).*(wx.*Rrr+vurr));

vgrf=cos(phi).*(vurf+wx.*Rrf)-wurf.*sin(phi);
vglf=cos(phi).*(vulf+wx.*Rlf)-wulf.*sin(phi);
vglr=cos(phi).*(vulr+wx.*Rlr)-wulr.*sin(phi);
vgrr=cos(phi).*(vurr+wx.*Rrr)-wurr.*sin(phi);
% vglr=Vslr(2)+wx*(lslr+Rlr);
% vgrr=Vsrr(2)+wx*(lsrr+Rrr);

% tire slip angle of each wheel
alpha_rf=atan(vgrf./ugrf)-delta;
alpha_lf=atan(vglf./uglf)-delta;
alpha_lr=atan(vglr./uglr);
alpha_rr=atan(vgrr./ugrr);

% % longitudinal slips
% if Rrf*wrf>(ugrf*cos(delta)+vgrf*sin(delta))
% % s_rf=(Rrf*wrf-(ugrf*cos(delta)+vgrf*sin(delta)))/(Rrf*wrf);
% kappa_rf=1-abs((ugrf*cos(delta)+vgrf*sin(delta))/(Rrf*wrf));
% else 
% % s_rf=(Rrf*wrf-(ugrf*cos(delta)+vgrf*sin(delta)))/(ugrf*cos(delta)+vgrf*sin(delta));
% kappa_rf=abs((Rrf*wrf)/(ugrf*cos(delta)+vgrf*sin(delta)))-1;
% end
% if Rlf*wlf>(uglf*cos(delta)+vglf*sin(delta))
% % s_lf=(Rlf*wlf-(uglf*cos(delta)+vglf*sin(delta)))/(Rlf*wlf);
% kappa_lf=1-abs((uglf*cos(delta)+vglf*sin(delta))/(Rlf*wlf));
% else
% % s_lf=(Rlf*wlf-(uglf*cos(delta)+vglf*sin(delta)))/(uglf*cos(delta)+vglf*sin(delta));    
% kappa_lf=abs((Rlf*wlf)/(uglf*cos(delta)+vglf*sin(delta)))-1;
% end
% if Rlr*wlr>uglr
% % s_lr=(Rlr*wlr-uglr)/(Rlr*wlr);
% kappa_lr=1-abs(uglr/(Rlr*wlr));
% else 
% % s_lr=(Rlr*wlr-uglr)/(uglr);
% kappa_lr=abs((Rlr*wlr)/uglr)-1;
% end
% if Rrr*wrr>ugrr
% % s_rr=(Rrr*wrr-ugrr)/(Rrr*wrr);
% kappa_rr=1-abs(ugrr/(Rrr*wrr));
% else
% % s_rr=(Rrr*wrr-ugrr)/(ugrr);
% kappa_rr=abs((Rrr*wrr)/ugrr)-1;
% end

kappa_lf =  -( 1 - (Rlf.*omega1./uglf) ) ;
kappa_rf =  -( 1 - (Rrf.*omega2./ugrf) ) ;
kappa_lr =  -( 1 - (Rlr.*omega3./uglr) ) ;
kappa_rr =  -( 1 - (Rrr.*omega4./ugrr) ) ;

% aerodynamic forces

dragForce = 0.5*airDens*dragCoeff*frontalArea*vx.^2;    
downforce       = 0.5*airDens*liftCoeff*frontalArea*vx.^2;
frontDownforce  = aeroBalance*downforce;
rearDownforce   = (1-aeroBalance)*downforce;

% get the tire loads once the tire planar forces depends on the loads
% to account for the wheel lift-off, when the tire radial compression becomes less than zero, the tire normal force Fz=0

Fzgrf=xtrf*ktf+frontDownforce/2;

Fzglf=xtlf*ktf+frontDownforce/2;

Fzglr=xtlr*ktr+rearDownforce/2;

Fzgrr=xtrr*ktr+rearDownforce/2;

% get forces and moments from the tire
% Fytlf
% Fytrf
% Fytlr
% Fytrr
% Fxtlf
% Fxtrf
% Fxtlr
% Fxtrr
% Mztlf
% Mztrf
% Mztlr
% Mztrr
[Fxtlf,Fytlf,Mztlf]  = Solver.MF5ss_eval(Fzglf,kappa_lf,alpha_lf,0,MF_f);    
[Fxtrf,Fytrf,Mztrf]  = Solver.MF5ss_eval(Fzgrf,kappa_rf,alpha_rf,0,MF_f);
[Fxtlr,Fytlr,Mztlr]  = Solver.MF5ss_eval(Fzglr,kappa_lr,alpha_lr,0,MF_r);    
[Fxtrr,Fytrr,Mztrr]  = Solver.MF5ss_eval(Fzgrr,kappa_rr,alpha_rr,0,MF_r);   

Fxglf=Fxtlf.*cos(delta)-Fytlf.*sin(delta);
Fxgrf=Fxtrf.*cos(delta)-Fytrf.*sin(delta);
Fxglr=Fxtlr;
Fxgrr=Fxtrr;
% Fxglr=Fxtlr*cos(delta)-Fytlr*sin(delta);
% Fxgrr=Fxtrr*cos(delta)-Fytrr*sin(delta);
Fyglf=Fxtlf.*sin(delta)+Fytlf.*cos(delta);
Fygrf=Fxtrf.*sin(delta)+Fytrf.*cos(delta);
% Fyglr=Fxtlr*sin(delta)+Fytlr*cos(delta);
% Fygrr=Fxtrr*sin(delta)+Fytrr*cos(delta);
Fyglr=Fytlr;
Fygrr=Fytrr;

% the tire force in coordinate frame 2
Fglf=[Fxglf;Fyglf;Fzglf];
Fgrf=[Fxgrf;Fygrf;Fzgrf];
Fglr=[Fxglr;Fyglr;Fzglr];
Fgrr=[Fxgrr;Fygrr;Fzgrr];

R_y = zeros(3,3,length(x(1,:))); R_x = zeros(3,3,length(x(1,:)));
Fgslf = zeros(3,length(x(1,:)));
Fgsrf = zeros(3,length(x(1,:)));
Fgslr = zeros(3,length(x(1,:)));
Fgsrr = zeros(3,length(x(1,:)));
%rotation matrix
for i=1:length(x(1,:))
    R_y(:,:,i) = [cos(theta(i)) 0 -sin(theta(i));0 1 0;sin(theta(i)) 0 cos(theta(i))];
    R_x(:,:,i) = [1 0 0;0 cos(phi(i)) sin(phi(i));0 -sin(phi(i)) cos(phi(i))];

    Fgslf(:,i)=R_x(:,:,i)*R_y(:,:,i)*Fglf(:,i);
    Fgsrf(:,i)=R_x(:,:,i)*R_y(:,:,i)*Fgrf(:,i);
    Fgslr(:,i)=R_x(:,:,i)*R_y(:,:,i)*Fglr(:,i);
    Fgsrr(:,i)=R_x(:,:,i)*R_y(:,:,i)*Fgrr(:,i);
end

% R_y=[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
% R_x=[1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];

% the force acting at the tire ground contact patch in coordinate 1 
% Fgslf=R_x*R_y*Fglf;
% Fgsrf=R_x*R_y*Fgrf;
% Fgslr=R_x*R_y*Fglr;
% Fgsrr=R_x*R_y*Fgrr;

% Fgslf=Fglf;
% Fgsrf=Fgrf;
% Fgslr=Fglr;
% Fgsrr=Fgrr;

%% the sprung mass dynamics module for each corner
% input: theta phi psi u v w wx wy wz(vehicle body), uuij, vuij, wuij, uuij_dot, vuij_dot, wuij_dot(unsprung mass)
% output: Fxsij Fysij Fzsij Mxij Myij Mzij(the sprung mass forces and moments at each corner)
% the forces transmitted to the sprung mass along x y z directions of coordinate 1
% right front corner
Fxsrf= Fgsrf(1,:)+muf*g*sin(theta)+muf*wz.*vurf-muf*wy.*wurf;
Fysrf= Fgsrf(2,:)-muf*g*sin(phi).*cos(theta)+muf*wx.*wurf-muf*wz.*uurf;
Fzsrf= xsrf*ksrf+dxsrf*bsrf;
% left front corner
Fxslf= Fgslf(1,:)+muf*g*sin(theta)+muf*wz.*vulf-muf*wy.*wulf;
Fyslf= Fgslf(2,:)-muf*g*sin(phi).*cos(theta)+muf*wx.*wulf-muf*wz.*uulf;
Fzslf= xslf*kslf+dxslf*bslf;
% left rear corner
Fxslr= Fgslr(1,:)+mur*g*sin(theta)+mur*wz.*vulr-mur*wy.*wulr;
Fyslr= Fgslr(2,:)-mur*g*sin(phi).*cos(theta)+mur*wx.*wulr-mur*wz.*uulr;
Fzslr= xslr*kslr+dxslr*bslr;
% right rear coner
Fxsrr= Fgsrr(1,:)+mur*g*sin(theta)+mur*wz.*vurr-mur*wy.*wurr;
Fysrr= Fgsrr(2,:)-mur*g*sin(phi).*cos(theta)+mur*wx.*wurr-mur*wz.*uurr;
Fzsrr= xsrr*ksrr+dxsrr*bsrr;
% the force represents the additional load transfer that occurs at the wheels
Fdzrf=(Fgsrf(2,:).*Rrf+Fysrf.*lsrf+Fgslf(2,:).*Rlf+Fyslf.*lslf-(Fysrf+Fyslf)*hrcf)/cf;
Fdzlf=-Fdzrf;
Fdzrr=(Fgsrr(2,:).*Rrr+Fysrr.*lsrr+Fgslr(2,:).*Rlr+Fyslr.*lslr-(Fysrr+Fyslr)*hrcr)/cr;
Fdzlr=-Fdzrr;
% the roll moment transmitted to the sprung mass
Mxrf=Fysrf*hrcf;
Mxlf=Fyslf*hrcf;
Mxlr=Fyslr*hrcr;
Mxrr=Fysrr*hrcr;
% the moments transmitted to the sprung mass by the suspension wy wz directions
Myrf=-(Fgsrf(1,:).*Rrf+Fxsrf.*lsrf);
Mylf=-(Fgslf(1,:).*Rlf+Fxslf.*lslf);
Mylr=-(Fgslr(1,:).*Rlr+Fxslr.*lslr);
Myrr=-(Fgsrr(1,:).*Rrr+Fxsrr.*lsrr);
Mzrf=Mztrf;
Mzlf=Mztlf;
Mzlr=Mztlr;
Mzrr=Mztrr;

% call powetrain calculation function
[Tfl,Tfr,Trl,Trr,~,~] = Solver.getWheelTorque14DOF(tau,p,Vehicle,x);

%% vehicle body dynamics module
% input: Fxsij Fysij Fzsij Mxij Myij Mzij(forces and moments transmitted to sprung mass at each corner)
% output: theta psi phi u v w wx wy wz u_dot v_dot w_dot wx_dot wy_dot wz_dot 
% u_dot=wz*v-wy*w+(1/m)*(Fxslf+Fxsrf+Fxslr+Fxsrr)+g*sin(theta);
u_dot=1/m*(Tfl+Tfr-Trl-Trr+Fxslf+Fxsrf+Fxslr+Fxsrr-dragForce+m*wz.*v-m*wy.*w+m*g*sin(theta));
v_dot=wx.*w-wz.*u+(1/m)*(Fyslf+Fysrf+Fyslr+Fysrr)-g*sin(phi).*cos(theta);
w_dot=wy.*u-wx.*v+(1/m)*(Fzslf+Fzsrf+Fzslr+Fzsrr+Fdzlf+Fdzrf+Fdzlr+Fdzrr)-g*cos(phi).*cos(theta);
wx_dot=(1/Jx)*((Mxlf+Mxrf+Mxlr+Mxrr)+(Fzslf-Fzsrf)*cf/2+(Fzslr-Fzsrr)*cr/2);
wy_dot=(1/Jy)*((Mylf+Myrf+Mylr+Myrr)+(Fzslr+Fzsrr)*b-(Fzslf+Fzsrf)*a);
wz_dot=(1/Jz)*((Mzlf+Mzrf+Mzlr+Mzrr)+(Fyslf+Fysrf)*a-(Fyslr+Fysrr)*b+(-Fxslf+Fxsrf)*cf/2 +(-Fxslr+Fxsrr)*cr/2);
% the cardan angles are obtained by performing the integration of the following quations
dtheta=wy.*cos(phi)-wz.*sin(phi);
dpsi=(wy.*sin(phi))./cos(theta)+(wz.*cos(phi))./cos(theta);
dphi=wx+wy.*sin(phi).*tan(theta)+wz.*cos(phi).*tan(theta);
% calculate the position in inertial coordinate
% Y_dot=u.*sin(psi)+v.*cos(psi);
% X_dot=u.*cos(psi)-v.*sin(psi);


% the unsprung mass vertical velocity wuij 
dzufl=(1/muf)*(cos(phi).*(cos(theta).*(Fzglf-muf*g)+sin(theta).*Fxglf)-sin(phi).*Fyglf-Fdzlf-xslf*kslf-dxslf*bslf-muf*(vulf.*wx-uulf.*wy));
dzufr=(1/muf)*(cos(phi).*(cos(theta).*(Fzgrf-muf*g)+sin(theta).*Fxgrf)-sin(phi).*Fygrf-Fdzrf-xsrf*ksrf-dxsrf*bsrf-muf*(vurf.*wx-uurf.*wy));
dzurl=(1/mur)*(cos(phi).*(cos(theta).*(Fzglr-mur*g)+sin(theta).*Fxglr)-sin(phi).*Fyglr-Fdzlr-xslr*kslr-dxslr*bslr-mur*(vulr.*wx-uulr.*wy));
dzurr=(1/mur)*(cos(phi).*(cos(theta).*(Fzgrr-mur*g)+sin(theta).*Fxgrr)-sin(phi).*Fygrr-Fdzrr-xsrr*ksrr-dxsrr*bsrr-mur*(vurr.*wx-uurr.*wy));



dx(1,:) = x(15,:);
dx(2,:) = x(16,:);
dx(3,:) = w;
dx(4,:) = dphi;
dx(5,:) = dtheta;
dx(6,:) = dpsi;
dx(7,:) = z1_dot;
dx(8,:) = z2_dot;
dx(9,:) = z3_dot;
dx(10,:) = z4_dot;
dx(11,:) = omega1;
dx(12,:) = omega2;
dx(13,:) = omega3;
dx(14,:) = omega4;
dx(15,:) = u_dot;
dx(16,:) = v_dot;
dx(17,:) = w_dot;
dx(18,:) = wx_dot;
dx(19,:) = wy_dot;
dx(20,:) = wz_dot;
dx(21,:) = dzufl;
dx(22,:) = dzufr;
dx(23,:) = dzurl;
dx(24,:) = dzurr;
dx(25,:) = (1/Iwheel_1)*(Tfl-Fxtlf*Rlf);
dx(26,:) = (1/Iwheel_2)*(Tfr-Fxtrf*Rrf);
dx(27,:) = (1/Iwheel_3)*(Trl-Fxtlr*Rlr);
dx(28,:) = (1/Iwheel_4)*(Trr-Fxtrr*Rrr);  
dx(29,:) = steerRate;
dx(30,:) = tauRate;

%% Equality Constraints : Tire Saturation
epsKap = 1e-5;
epsAlp = 1e-5;

[Fx1_t,~,~]     = Solver.MF5ss_eval(Fzglf,kappa_lf+epsKap,alpha_lf,0,MF_f);    
[Fx2_t,~,~]     = Solver.MF5ss_eval(Fzgrf,kappa_rf+epsKap,alpha_rf,0,MF_f);       
[Fx3_t,~,~]     = Solver.MF5ss_eval(Fzglr,kappa_lr+epsKap,alpha_lr,0,MF_r);    
[Fx4_t,~,~]     = Solver.MF5ss_eval(Fzgrr,kappa_rr+epsKap,alpha_rr,0,MF_r);    

[~,Fy1_t,~]     = Solver.MF5ss_eval(Fzglf,kappa_lf,alpha_lf+epsKap,0,MF_f);    
[~,Fy2_t,~]     = Solver.MF5ss_eval(Fzgrf,kappa_rf,alpha_rf+epsKap,0,MF_f);
[~,Fy3_t,~]     = Solver.MF5ss_eval(Fzglr,kappa_lr,alpha_lr+epsKap,0,MF_r);    
[~,Fy4_t,~]     = Solver.MF5ss_eval(Fzgrr,kappa_rr,alpha_rr+epsKap,0,MF_r);   


% Calculation of Slip Stiffness
C_x_1 = (Fx1_t - Fx1)/epsKap;
C_x_2 = (Fx2_t - Fx2)/epsKap;
C_x_3 = (Fx3_t - Fx3)/epsKap;
C_x_4 = (Fx4_t - Fx4)/epsKap;

% Calculation of Cornering Stiffness
C_y_1 = (Fy1_t - Fy1)/epsAlp;
C_y_2 = (Fy2_t - Fy2)/epsAlp;
C_y_3 = (Fy3_t - Fy3)/epsAlp;
C_y_4 = (Fy4_t - Fy4)/epsAlp;   

% Groups the inequality constraint into a single vector
%     g_ineq = [C_x_1;
%      C_x_2;
%      C_y_1;
%      C_y_2;
%      C_y_3;
%      C_y_4;
%      C_x_3;
%      C_x_4];

%% Equality Constraints : Vehicle Stability

% Understeer angle
wb = a+b;
delta_kin = wb*wz./vx;
theta_uang = delta-delta_kin;
theta_lim_uang = deg2rad(8);
C_uang = (theta_uang/theta_lim_uang).^2-1;

% Bounding the maximum value of beta
beta = atan(v./u);
beta_lim = deg2rad(5);
C_beta = (beta/beta_lim).^2-1;

%     yaw_stiffness = jacobian(YMnet,beta);   %yaw stiffness
%     gradient(dot(A,A),A)
%         yaw_stiffness = evalf(gradient(YMnet))./evalf(gradient(beta));   %yaw stiffness

%     g_ineq = -[C_x_1;
%      C_x_2;
%      C_x_3;
%      C_x_4];

g_ineq = -[C_x_1;
 C_x_2;
 C_x_3;
 C_x_4;...
 -C_uang;...
 -C_beta];


%% Saving Constraints : Fuel Usage

mass_flow  = Solver.mass_flow_calculation(x,tau,Vehicle);
fuel_usage = ds*(Sf.*mass_flow)*ones(length(u),1) ;

%% Saving Constraints : Tire Energy
  
% Longitudinal Sliding Speed
long_sliding_speed_1 = uglf.*kappa_lf;
long_sliding_speed_2 = ugrf.*kappa_rf;
long_sliding_speed_3 = uglr.*kappa_lr;
long_sliding_speed_4 = ugrr.*kappa_rr;

% Lateral Sliding Speed
lat_sliding_speed_1 = uglf.*tan(alpha_lf);
lat_sliding_speed_2 = ugrf.*tan(alpha_rf);
lat_sliding_speed_3 = uglr.*tan(alpha_rf);
lat_sliding_speed_4 = ugrr.*tan(alpha_rr);

% IMPORTANT
% The way the tire energies are being calculated are not entirely
% correct now, because I'm taking the integral over the tire power over
% the distance instead of time. Because the speed is one of the states
% of the problem and therefore it is a casadi MX variable, I cannot
% calculate the (variable) time step

sliding_power_lateral_1 = lat_sliding_speed_1 .* Fytlf ; %[w]
%     sliding_energy_lateral_1 = cumtrapz(cumLaptime,abs(sliding_power_lateral_1))./1; % J
sliding_energy_lateral_1 = ds*abs(sliding_power_lateral_1)*ones(length(n),1); % J

sliding_power_lateral_2 = lat_sliding_speed_2 .* Fytrf ; %[w]
%     sliding_energy_lateral_2 = cumtrapz(cumLaptime,abs(sliding_power_lateral_2))./1; % J
sliding_energy_lateral_2 = ds*abs(sliding_power_lateral_2)*ones(length(n),1); % J

sliding_power_lateral_3 = lat_sliding_speed_3 .* Fytlr ; %[w]
%     sliding_energy_lateral_3 = cumtrapz(cumLaptime,abs(sliding_power_lateral_3))./1; % J
sliding_energy_lateral_3 = ds*abs(sliding_power_lateral_3)*ones(length(n),1); % J

sliding_power_lateral_4 = lat_sliding_speed_4 .* Fytrr ; %[w]
%     sliding_energy_lateral_4 = cumtrapz(cumLaptime,abs(sliding_power_lateral_4))./1; % J 
sliding_energy_lateral_4 = ds*abs(sliding_power_lateral_4)*ones(length(n),1); % J 

sliding_power_longitudinal_1 = long_sliding_speed_1 .* Fxtlf ; %[w]
sliding_energy_longitudinal_1 = ds*abs(sliding_power_longitudinal_1)*ones(length(n),1); % J

sliding_power_longitudinal_2 = long_sliding_speed_2 .* Fxtrf ; %[w]
sliding_energy_longitudinal_2 = ds*abs(sliding_power_longitudinal_2)*ones(length(n),1); % J

sliding_power_longitudinal_3 = long_sliding_speed_3 .* Fxtlr ; %[w]
sliding_energy_longitudinal_3 = ds*abs(sliding_power_longitudinal_3)*ones(length(n),1); % J

sliding_power_longitudinal_4 = long_sliding_speed_4 .* Fxtrr ; %[w]
sliding_energy_longitudinal_4 = ds*abs(sliding_power_longitudinal_4)*ones(length(n),1); % J

% Combined Energy Dissipated into the 4 Tires

fl_tire_energy_combined = sliding_energy_lateral_1 + ...
    sliding_energy_longitudinal_1;

fr_tire_energy_combined = sliding_energy_lateral_2 + ...
    sliding_energy_longitudinal_2;

rl_tire_energy_combined = sliding_energy_lateral_3 + ...
    sliding_energy_longitudinal_3;

rr_tire_energy_combined = sliding_energy_lateral_4 + ...
    sliding_energy_longitudinal_4;    



%% saving constraints  


    saving_constraints.tire_energy = fl_tire_energy_combined + fr_tire_energy_combined +...
        rl_tire_energy_combined + rr_tire_energy_combined;

% saving_constraints.tire_energy = fl_tire_energy_combined;

saving_constraints.fuel = fuel_usage(end);