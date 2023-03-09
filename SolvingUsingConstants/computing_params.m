clear, clc

%% Compute Ksh

data_speed_step = load('Data/SpeedSteps_Noise_U.mat');
AUVparams = loadjson('Conf/AUVParameters.json');

u = data_speed_step.vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;
Fx = [data_speed_step.forces_values.HydrodynamicForces.Fx_N.Data(:)];

% figure(1)
% subplot(1, 2, 1)
% plot(Fx)
% title('u')
% subplot(1, 2, 2)
% plot(u)
% title('Fx')

rho = AUVparams.Environment.Rho;
nu = AUVparams.Environment.nu;
Sref = AUVparams.Mecanic.surface_reference;
Lref = AUVparams.Mecanic.length_reference;

U = 1/2 * rho * Sref * 0.075 .* u .* abs(u) ./ (log10(u * Lref / nu) - 2).^2;

%Ksh = U \ Fx;  
Ksh = lsqr(U, Fx);

AUVparams.Hydrodynamic.Ksh = Ksh;

%% Compute CZuw and CZuq

data_pitch_step = load('Data/PitchSteps_NoiseOn_Pitch_W_Q.mat');

u = data_pitch_step.vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;
w = data_pitch_step.vehicle_state.nu.rAUV_WaterSpeed.Wwater_ms.Data;
q = data_pitch_step.vehicle_state.nu.AngularSpeed.Q_rads.Data;

CZ0 = AUVparams.Hydrodynamic.CZ0;
CM0 = AUVparams.Hydrodynamic.CM0;

Fz = [data_pitch_step.forces_values.HydrodynamicForces.Fz_N.Data(:)];
My = [data_pitch_step.forces_values.HydrodynamicForces.My_Nm.Data(:)];

F = Fz - 1/2 * rho * Sref * CZ0 .* u .* abs(u) ;
A = 1/2 * rho * Sref .* u .* w;
B = 1/2 * rho * Sref * Lref .* u .* q;

X = [A B]\F;
AUVparams.Hydrodynamic.CZuw = X(1);
AUVparams.Hydrodynamic.CZuq = X(2);

M = My - 1/2 * rho * Sref * Lref * CM0 .* u .* abs(u) ;
A = 1/2 * rho * Sref * Lref .* u .* w;
B = 1/2 * rho * Sref * Lref^2 .* u .* q;

X = [A B]\M;

AUVparams.Hydrodynamic.CMuw = X(1);
AUVparams.Hydrodynamic.CMuq = X(2);

%% Compute CYuv, CYur, CNuv and CNur

data_dieudonne = load('Data/DieuDonne_AStep_NoiseOn_U_V_R_A.mat');

u = data_dieudonne.vehicle_state.nu.rAUV_WaterSpeed.Uwater_ms.Data;
v = data_dieudonne.vehicle_state.nu.rAUV_WaterSpeed.Vwater_ms.Data;
r = data_dieudonne.vehicle_state.nu.AngularSpeed.R_rads.Data;

Fy = [data_dieudonne.forces_values.HydrodynamicForces.Fy_N.Data(:)];
Mz = [data_dieudonne.forces_values.HydrodynamicForces.Mz_Nm.Data(:)];

A = 1/2 * rho * Sref .* u .* v;
B = 1/2 * rho * Sref * Lref .* u .* r;

% Define the function
% fxn = @(C) (Fy - (A .* C(1) + B .* C(2)) ).^2;

% Defining the constraint for the fmincon
CZuw = AUVparams.Hydrodynamic.CZuw;
CZuq = AUVparams.Hydrodynamic.CZuq;

x0 = [CZuw; -CZuq]; % Initial values

UB = [0.8*CZuw; -1.025*CZuq];
LB = [1.2*CZuw; -0.975*CZuq];

X = fminsearchcon(@(x) sum((Fy - (A .* x(1) + B .* x(2)) ).^2), x0, LB, UB);
AUVparams.Hydrodynamic.CYuv = X(1);
AUVparams.Hydrodynamic.CYur = X(2);

% Computing CNuv and CNur
A = 1/2 * rho * Sref * Lref .* u .* v;
B = 1/2 * rho * Sref * Lref^2 .* u .* r;

CMuw = AUVparams.Hydrodynamic.CMuw;
CMuq = AUVparams.Hydrodynamic.CMuq;

x0 = [-CMuw; CMuq]; % Initial values

UB = [-1.1*CMuw; 0.9*CMuq];
LB = [-0.9*CMuw; 1.1*CMuq];

X = fminsearchcon(@(x) sum((Mz - (A .* x(1) + B .* x(2)) ).^2), x0, LB, UB);
AUVparams.Hydrodynamic.CNuv = X(1);
AUVparams.Hydrodynamic.CNur = X(2);

%% Exporting the Params

savejson('', AUVparams, 'Conf/AUVParameters.json');
