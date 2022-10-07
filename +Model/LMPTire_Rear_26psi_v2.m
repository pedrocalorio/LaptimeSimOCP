function [mfCoeff] = LMPTire_Rear_26psi_v2(~)
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
    MF.scaling.lmux = 0.9;   % longitudinal friction scaling
    MF.scaling.lKx  = 0.9;   % longitudinal slip stiffness scaling
    MF.scaling.lmuy = 0.85;   % lateral friction scaling
    MF.scaling.lKy  = 0.9;   % cornering stiffness scaling
    MF.scaling.lgay = 1;   % scaling of camber influence on lateral force
    MF.scaling.ltr  = 1;   % scaling of pneumatic trail
    MF.scaling.lres = 1;   % scaling of residual moment 
    MF.scaling.lgaz = 1;   % scaling of camber influence on self aligning moment
    MF.scaling.ls   = 1;   % scaling of moment arm Fx (combined slip)

    % longitudinal force coefficients
    MF.long.pCx1 =  1.57;
    MF.long.pDx1 =  2.06;
    MF.long.pDx2 = -0.149;
    MF.long.pEx1 = -0.795;
    MF.long.pEx2 =  -0.314;
    MF.long.pEx3 =  1.51;
    MF.long.pKx1 = 48.7;
    MF.long.pKx2 =  -0.906;
    MF.long.pKx3 =  0.186;
    MF.long.rBx1 = 13.6;
    MF.long.rBx2 = 12.5;
    MF.long.rCx1 =  0.941;

    % lateral force coefficients
    MF.lat.pCy1 =  1.46;
    MF.lat.pDy1 =  1.81*0.85;
%     MF.lat.pDy2 = -0.619;
    MF.lat.pDy2 = -0.619*1.0;
    MF.lat.pDy3 =  1.43;
    MF.lat.pEy1 = -0.739;
    MF.lat.pEy2 = -1.53;
    MF.lat.pEy4 = -7.28;
    MF.lat.pKy1 = 38.3;
    MF.lat.pKy2 =  1.99;
    MF.lat.pKy3 =  -0.142;
    MF.lat.pHy3 =  0.0105;
    MF.lat.pVy3 =  2.22;
    MF.lat.pVy4 =  2.2;
    MF.lat.rBy1 =  15.1;
    MF.lat.rBy2 =  11.2;
    MF.lat.rCy1 =  1.01;

    % self aligning moment coefficients
    MF.align.qBz1  =  5.57;
    MF.align.qBz2  = 4.15;
    MF.align.qBz3  = -0.77;
    MF.align.qBz5  =  -0.736;
    MF.align.qBz9  =  28.9;
    MF.align.qBz10 =  -0.17;
    MF.align.qCz1  =  1.08;
    MF.align.qDz1  =  0.13;
    MF.align.qDz2  = 0.0221;
    MF.align.qDz4  = -38.5;
    MF.align.qDz8  =  1.72;
    MF.align.qDz9  =  1.01;
    MF.align.qEz1  = -21.6;
    MF.align.qEz2  = 71.8;
    MF.align.qEz3  =  -118;
    MF.align.qEz5  = -16.7;
    MF.align.qHz3  =  0.262;
    MF.align.qHz4  = -0.00356;
    MF.align.ssz2  =  0.0384;
    MF.align.ssz3  = -1.73;
    MF.align.ssz4  =  1.64;


    mfCoeff = MF;
end

