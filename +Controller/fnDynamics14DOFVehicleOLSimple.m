function dx = fnDynamics14DOFVehicleOLSimple(t,x,u_c,p,Vehicle)

addpath('C:/dev/libraries/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

% kappa = interp1(Track.distance,Track.curv,Track.sLap,'spline');

% Distance Step
% ds = (Track.sLap(end) - Track.sLap(1)) / (length(x(1,:)) - 1);

% Gravity constant
g = 9.81;

%% System states


% x = x(1); % longitudinal displacement
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
tauRate   = u_c(2,:);

%% Calculating the cumulative laptime

% Sf = (1 - n.*kappa)./(u.*cos(xi)-v.*sin(xi));
% cumLaptime = ds*Sf*ones(length(n),1);

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
% xtirf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
% xtilf=((m*g*b)/(2*(a+b))+muf*g)/ktf;
% xtilr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
% xtirr=((m*g*a)/(2*(a+b))+mur*g)/ktr;
xtlf = -z1;
xtrf = -z2;
xtlr = -z3;
xtrr = -z4;

%the initial spring compression xsif 
% xsirf=(m*g*b)/(2*(a+b)*ksrf);
% xsilf=(m*g*b)/(2*(a+b)*kslf);
% xsilr=(m*g*a)/(2*(a+b)*kslr);
% xsirr=(m*g*a)/(2*(a+b)*ksrr);
xslf=z1 + (heave + a*sin(theta) + cf/2*sin(phi) );
xsrf=z2 + (heave + a*sin(theta) - cf/2*sin(phi) );
xslr=z3 + (heave - b*sin(theta) + cr/2*sin(phi) );
xsrr=z4 + (heave - b*sin(theta) - cr/2*sin(phi) );

% the initial length of the strut
% lsirf=h-(r0-xtirf);
% lsilf=h-(r0-xtilf);
% lsilr=h-(r0-xtilr);
% lsirr=h-(r0-xtirr);
lslf=h-(r01+xtlf+xslf);
lsrf=h-(r02+xtrf+xsrf);
lslr=h-(r03+xtlr+xslr);
lsrr=h-(r04+xtrr+xsrr);

% the transforming matrix 
Mrf=[0 0 cf/2;0 0 a;-cf/2 -a 0];
Mlf=[0 0 -cf/2;0 0 a;cf/2 -a 0];
Mlr=[0 0 -cr/2;0 0 -b;cr/2 b 0];
Mrr=[0 0 cr/2;0 0 -b;-cr/2 b 0];

% initial velocity of strut mounting point in x y z by transforming the c.g. velocities
wm=[wx;wy;wz];
vm=[u;v;w];

% Vsrf = zeros(3,length(n)); Vslf = zeros(3,length(n));
% Vsrr = zeros(3,length(n)); Vslr = zeros(3,length(n));
% Vsrf = []; Vslf = [];
% Vsrr = []; Vslr = [];
% % vou precisar fazer um for loop
% for i=1:length(x(1,:))
%     Vsrf_i=Mrf*wm(:,i)+vm(:,i); % front right mounting point x y z velocity in coordinate frame 1
%     Vslf_i=Mlf*wm(:,i)+vm(:,i); % front left
%     Vslr_i=Mlr*wm(:,i)+vm(:,i); % rear left
%     Vsrr_i=Mrr*wm(:,i)+vm(:,i); % rear right
% 
%     Vsrf=[Vsrf Vsrf_i];
%     Vslf=[Vslf Vslf_i];
%     Vslr=[Vslr Vslr_i];
%     Vsrr=[Vsrr Vsrr_i];
% end
Vslf=Mlf*wm+vm; % front left
Vsrf=Mrf*wm+vm; % front right mounting point x y z velocity in coordinate frame 1
Vslr=Mlr*wm+vm; % rear left
Vsrr=Mrr*wm+vm; % rear right

% initial unsprung mass velocity
uulf=Vslf(1,:)-lslf.*wy; 
uurf=Vsrf(1,:)-lsrf.*wy; 
uulr=Vslr(1,:)-lslr.*wy; 
uurr=Vsrr(1,:)-lsrr.*wy; 
vulf=Vslf(2,:)+lslf.*wx; 
vurf=Vsrf(2,:)+lsrf.*wx; 
vulr=Vslr(2,:)+lslr.*wx; 
vurr=Vsrr(2,:)+lsrr.*wx;

wulf=z1_dot;
wurf=z2_dot;
wulr=z3_dot;
wurr=z4_dot;

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

Rlf=(r01-xtlf)./(cos(theta).*cos(phi));

Rrf=(r02-xtrf)./(cos(theta).*cos(phi));

Rlr=(r03-xtlr)./(cos(theta).*cos(phi));

Rrr=(r04-xtrr)./(cos(theta).*cos(phi));

% the longitudinal and lateral velocities at the tire contact patch in coordinate frame 2
uglf=cos(theta).*(uulf-wy.*Rlf)+sin(theta).*(wulf.*cos(phi)+sin(phi).*(wx.*Rlf+vulf));
ugrf=cos(theta).*(uurf-wy.*Rrf)+sin(theta).*(wurf.*cos(phi)+sin(phi).*(wx.*Rrf+vurf));
uglr=cos(theta).*(uulr-wy.*Rlr)+sin(theta).*(wulr.*cos(phi)+sin(phi).*(wx.*Rlr+vulr));
ugrr=cos(theta).*(uurr-wy.*Rrr)+sin(theta).*(wurr.*cos(phi)+sin(phi).*(wx.*Rrr+vurr));

vglf=cos(phi).*(vulf+wx.*Rlf)-wulf.*sin(phi);
vgrf=cos(phi).*(vurf+wx.*Rrf)-wurf.*sin(phi);
vglr=cos(phi).*(vulr+wx.*Rlr)-wulr.*sin(phi);
vgrr=cos(phi).*(vurr+wx.*Rrr)-wurr.*sin(phi);
% vglr=Vslr(2)+wx*(lslr+Rlr);
% vgrr=Vsrr(2)+wx*(lsrr+Rrr);

% tire slip angle of each wheel
alpha_lf=atan(vglf./uglf)-delta;
alpha_rf=atan(vgrf./ugrf)-delta;
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

% kappa_lf =  -( 1 - (Rlf.*omega1./uglf) ) ;
% kappa_rf =  -( 1 - (Rrf.*omega2./ugrf) ) ;
% kappa_lr =  -( 1 - (Rlr.*omega3./uglr) ) ;
% kappa_rr =  -( 1 - (Rrr.*omega4./ugrr) ) ;

kappa_lf =  (Rlf*omega1-(uglf*cos(delta)+vglf*sin(delta)))/(uglf*cos(delta)+vglf*sin(delta)) ;
kappa_rf =  (Rrf*omega2-(ugrf*cos(delta)+vgrf*sin(delta)))/(ugrf*cos(delta)+vgrf*sin(delta)) ;
kappa_lr =  (Rlr*omega3-uglr)/(uglr) ;
kappa_rr =  (Rrr*omega4-ugrr)/(ugrr) ;

% aerodynamic forces

dragForce       = 0.5*airDens*dragCoeff*frontalArea*u.^2;    
downforce       = 0.5*airDens*liftCoeff*frontalArea*u.^2;
frontDownforce  = aeroBalance*downforce;
rearDownforce   = (1-aeroBalance)*downforce;

% get the tire loads once the tire planar forces depends on the loads
% to account for the wheel lift-off, when the tire radial compression becomes less than zero, the tire normal force Fz=0

Fzglf=xtlf*ktf+frontDownforce/2;

Fzgrf=xtrf*ktf+frontDownforce/2;

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
[Fxtlf,Fytlf,Mztlf]  = Solver.MF5ss_eval(Fzglf,kappa_lf,-alpha_lf,0,MF_f);    
[Fxtrf,Fytrf,Mztrf]  = Solver.MF5ss_eval(Fzgrf,kappa_rf,-alpha_rf,0,MF_f);
[Fxtlr,Fytlr,Mztlr]  = Solver.MF5ss_eval(Fzglr,kappa_lr,-alpha_lr,0,MF_r);    
[Fxtrr,Fytrr,Mztrr]  = Solver.MF5ss_eval(Fzgrr,kappa_rr,-alpha_rr,0,MF_r);   

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

% R_y = zeros(3,3,length(x(1,:))); R_x = zeros(3,3,length(x(1,:)));
% Fgslf = zeros(3,length(x(1,:)));
% Fgsrf = zeros(3,length(x(1,:)));
% Fgslr = zeros(3,length(x(1,:)));
% Fgsrr = zeros(3,length(x(1,:)));

% Fgslf = [];
% Fgsrf = [];
% Fgslr = [];
% Fgsrr = [];
% %rotation matrix
% for i=1:length(x(1,:))
%     R_y = [cos(theta(i)) 0 -sin(theta(i));0 1 0;sin(theta(i)) 0 cos(theta(i))];
%     R_x = [1 0 0;0 cos(phi(i)) sin(phi(i));0 -sin(phi(i)) cos(phi(i))];
% 
%     Fgslf_i=R_x*R_y*Fglf(:,i);
%     Fgsrf_i=R_x*R_y*Fgrf(:,i);
%     Fgslr_i=R_x*R_y*Fglr(:,i);
%     Fgsrr_i=R_x*R_y*Fgrr(:,i);
% 
%     Fgslf = [Fgslf Fgslf_i];
%     Fgsrf = [Fgsrf Fgsrf_i];
%     Fgslr = [Fgslr Fgslr_i];
%     Fgsrr = [Fgsrr Fgsrr_i];
% end

R_y=[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
R_x=[1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];

%the force acting at the tire ground contact patch in coordinate 1 
Fgslf=R_x*R_y*Fglf;
Fgsrf=R_x*R_y*Fgrf;
Fgslr=R_x*R_y*Fglr;
Fgsrr=R_x*R_y*Fgrr;



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
[Tfl,Tfr,Trl,Trr,~,~] = Solver.getWheelTorque2(tau,p,Vehicle,omega3,omega4);

%% vehicle body dynamics module
% input: Fxsij Fysij Fzsij Mxij Myij Mzij(forces and moments transmitted to sprung mass at each corner)
% output: theta psi phi u v w wx wy wz u_dot v_dot w_dot wx_dot wy_dot wz_dot 
% u_dot=wz*v-wy*w+(1/m)*(Fxslf+Fxsrf+Fxslr+Fxsrr)+g*sin(theta);
u_dot=1/m*(Fxslf+Fxsrf+Fxslr+Fxsrr-dragForce+m*wz.*v-m*wy.*w+m*g*sin(theta));
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


dx = [u;
      v;
      w;
      dphi;
      dtheta;
      dpsi;
      z1_dot;
      z2_dot;
      z3_dot;
      z4_dot;
      omega1;
      omega2;
      omega3;
      omega4;
      u_dot;
      v_dot;
      w_dot;
      wx_dot;
      wy_dot;
      wz_dot;
      dzufl;
      dzufr;
      dzurl;
      dzurr;
      (1/Iwheel_1)*(Tfl-Fxglf.*Rlf);
      (1/Iwheel_2)*(Tfr-Fxgrf.*Rlr);
      (1/Iwheel_3)*(Trl-Fxglr.*Rrf);
      (1/Iwheel_4)*(Trr-Fxgrr.*Rrr);
      steerRate;
      tauRate];

