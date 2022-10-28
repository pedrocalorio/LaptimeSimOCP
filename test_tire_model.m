problem = struct();

problem = Initialize.fnInitVehicle(problem);

% range of alphas

alpha_range = linspace(-15,15,100);
alpha_range = deg2rad(alpha_range);

fz_range = linspace(1e3,10e3,100);

% range of kappas

kappa_range = linspace(-0.9,0.9,100);

Fz = 6e3;
gamma = 0;

MF_f = problem.dsSystem.vd.tire_1.MF;

MF_r = problem.dsSystem.vd.tire_3.MF;

for i=1:length(alpha_range)
    alpha = alpha_range(i);
    for j=1:length(kappa_range)
        kappa = kappa_range(j);
        [Fx_f(i,j),Fy_f(i,j),~] = Solver.MF5ss_eval(Fz,kappa,alpha,gamma,MF_f);
        [Fx_r(i,j),Fy_r(i,j),~] = Solver.MF5ss_eval(Fz,kappa,alpha,gamma,MF_r);
    end
end

for i=1:length(alpha_range)
    alpha = alpha_range(i);
        [Fx_fP(i),Fy_fP(i),~] = Solver.MF5ss_eval(Fz,0,alpha,gamma,MF_f);
        [Fx_rP(i),Fy_rP(i),~] = Solver.MF5ss_eval(Fz,0,alpha,gamma,MF_r);
end

for i=1:length(fz_range)
    Fz = fz_range(i);
        [~,Fy_fP_fz(i),~] = Solver.MF5ss_eval(Fz,0,deg2rad(5),gamma,MF_f);
        [~,Fy_rP_fz(i),~] = Solver.MF5ss_eval(Fz,0,deg2rad(5),gamma,MF_r);
end

for i=1:length(fz_range)
    Fz = fz_range(i);
        [Fx_fP_fz(i),~,~] = Solver.MF5ss_eval(Fz,0.05,deg2rad(0),gamma,MF_f);
        [Fx_rP_fz(i),~,~] = Solver.MF5ss_eval(Fz,0.05,deg2rad(0),gamma,MF_r);
end


surf(kappa_range,rad2deg(alpha_range),Fy_f)
hold on
surf(kappa_range,rad2deg(alpha_range),Fy_r)
ylabel('slip angle [deg]')
xlabel('slip ratio [-]')
zlabel('Fy')

figure;
plot(rad2deg(alpha_range),Fy_fP)
hold on
plot(rad2deg(alpha_range),Fy_rP)
legend('front','rear')
title('Fy vs SA')

%%
figure;
plot((fz_range),Fy_fP_fz/Fz)
hold on
plot((fz_range),Fy_rP_fz/Fz)
legend('front','rear')
title('\mi_{y} vs Fz','Interpreter','latex')

figure;
plot((fz_range),Fx_fP_fz/Fz)
hold on
plot((fz_range),Fx_rP_fz/Fz)
legend('front','rear')
title('\mi_{x} vs Fz','Interpreter','latex')

% figure
% surf(kappa_range,rad2deg(alpha_range),Fx_f)
% ylabel('slip angle [deg]')
% xlabel('slip ratio [-]')
% zlabel('Fx')