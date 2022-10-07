function [mfCoeff] = LMPTire_Front_26psi(~)
    % Tyre model Magic Formula coefficients
    %
    % Note 1: these coefficients use the Pacejka sign convention (See book: Tire and Vehicle Dynamics)
    % They are not compatible with the ISO sign convention used by TNO or ADAMS tyre models!
    % see tyre sign convention, Vehicle Dynamics lecture notes (4L150), page 118.
    %
    % Note 2: simplified, symmetric Magic Formula expressions are used, so a  number of coefficients 
    % have been eliminated (e.g. coefficients controlling static offsets: pHx1, pHx2, pHy1, pHy2, pVy1, pVy2, etc.)

    MF.fittyp       = 5;         % Magic Formula version [-]
    MF.R0           = 0.3300;    % unloaded tyre radius used in Magic Formula to make coefficients dimensionless [m]
    MF.Fz0          = 6000;      % reference force value (typically corresponding to the static vertical force), 
                                        % used in Magic Formula to make coefficients dimensionless [N]

    % limitations on the inputs of the Magic Formula (to prevent extrapolation errors)
    MF.limits.Fz    = [200 11564];
    MF.limits.kappa = [-1.5000 1.5000];
    MF.limits.alpha = [-1.5708 1.5708];
    MF.limits.gamma = [-0.2610 0.2610];

    % dimensionless scaling coefficients to modify global tyre behaviour
    MF.scaling.lmux = 0.8;   % longitudinal friction scaling
    MF.scaling.lKx  = 0.8;   % longitudinal slip stiffness scaling
    MF.scaling.lmuy = 0.8;   % lateral friction scaling
    MF.scaling.lKy  = 0.8;   % cornering stiffness scaling
    MF.scaling.lgay = 1;   % scaling of camber influence on lateral force
    MF.scaling.ltr  = 1;   % scaling of pneumatic trail
    MF.scaling.lres = 1;   % scaling of residual moment 
    MF.scaling.lgaz = 1;   % scaling of camber influence on self aligning moment
    MF.scaling.ls   = 1;   % scaling of moment arm Fx (combined slip)

    % longitudinal force coefficients
    MF.long.pCx1 =  1.55;
    MF.long.pDx1 =  1.99; 
    MF.long.pDx2 = 0.16;
    MF.long.pEx1 = -0.98;
    MF.long.pEx2 =  0.0275;
    MF.long.pEx3 =  1.91;
    MF.long.pKx1 = 50.9;
    MF.long.pKx2 =  -17.2;
    MF.long.pKx3 =  0.467;
    MF.long.rBx1 = 15.6;
    MF.long.rBx2 = 14.6;
    MF.long.rCx1 =  0.808;

    % lateral force coefficients
    MF.lat.pCy1 =  1.61;
    MF.lat.pDy1 =  1.56;
    MF.lat.pDy2 = -0.372;
    MF.lat.pDy3 =  1.53;
    MF.lat.pEy1 = 0.00104;
    MF.lat.pEy2 = 0.013;
    MF.lat.pEy4 = 79.9;
    MF.lat.pKy1 = 50.4;
    MF.lat.pKy2 =  3.13;
    MF.lat.pKy3 =  0.26;
    MF.lat.pHy3 =  0.0184;
    MF.lat.pVy3 =  1.77;
    MF.lat.pVy4 =  2.7;
    MF.lat.rBy1 =  15.1;
    MF.lat.rBy2 =  7.89;
    MF.lat.rCy1 =  1.00;

    % self aligning moment coefficients
    MF.align.qBz1  =  7.18;
    MF.align.qBz2  = 5.18;
    MF.align.qBz3  = -0.559;
    MF.align.qBz5  =  -2.87;
    MF.align.qBz9  =  5.48;
    MF.align.qBz10 =  0.297;
    MF.align.qCz1  =  1.18;
    MF.align.qDz1  =  0.111;
    MF.align.qDz2  = 0.0161;
    MF.align.qDz4  = -23;
    MF.align.qDz8  =  1.33;
    MF.align.qDz9  =  1.54;
    MF.align.qEz1  = -8.85;
    MF.align.qEz2  = 34.7;
    MF.align.qEz3  =  -56.2;
    MF.align.qEz5  = -11.3;
    MF.align.qHz3  =  0.114;
    MF.align.qHz4  = -0.42;
    MF.align.ssz2  =  0.0294;
    MF.align.ssz3  = -1.67;
    MF.align.ssz4  =  1.31;


    mfCoeff = MF;
end

