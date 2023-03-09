clear; clc; close all

%% Inputs

u_val = 10;
z_val = -50;

%% Computations

AUVparams = loadjson('Conf/AUVParameters.json');

rho = AUVparams.Environment.Rho;
nu = AUVparams.Environment.nu;
Sref = AUVparams.Mecanic.surface_reference;
Lref = AUVparams.Mecanic.length_reference;

CMuw = AUVparams.Hydrodynamic.CMuw;
CMuq = AUVparams.Hydrodynamic.CMuq;

CZuw = AUVparams.Hydrodynamic.CZuw;
CZuq = AUVparams.Hydrodynamic.CZuq;

CZ0 = AUVparams.Hydrodynamic.CZ0;
CM0 = AUVparams.Hydrodynamic.CM0;

m = AUVparams.Mecanic.mass;
Zw_dot = AUVparams.Mecanic.added_water_matrix(3, 3);
Zq_dot = AUVparams.Mecanic.added_water_matrix(3, 5);
Mv_dot = AUVparams.Mecanic.added_water_matrix(5, 2);
Iy = AUVparams.Mecanic.inertia_matrix(2, 2);
Mq_dot = AUVparams.Mecanic.added_water_matrix(5, 5);
zg = AUVparams.Mecanic.center_of_gravity_veh_ref(3);
zb = AUVparams.Mecanic.center_of_buoyancy_veh_ref(3);

CZ = AUVparams.XRearHelms.CZ;
CM = AUVparams.XRearHelms.CM;

gravity = 9.81;
V = AUVparams.Mecanic.volume;

W = m * gravity;
B = V*rho*gravity;

k = 0.5*rho*Sref;

syms u0 w q z int_z theta BAR

A_cal = u0 * [0; m*w];
B_cal = [0, m*(zg*q + u0); -m*(zg*q + u0), 0] * [w; q];
C_cal = k * [CZ0*abs(u0); Lref * CM0 * abs(u0)] * u0;
D_cal = k * [CZuw * u0, CZuq*Lref*u0; CMuw*Lref*u0, Lref^2*CMuq*u0]*[w; q];
E_cal = [(W - B)*cos(theta); -(zg*W - zb*B)*sin(theta)];
F_cal = k * [CZ*BAR*u0* abs(u0); Lref*CM*BAR*u0*abs(u0)];
M_cal = [(m+Zw_dot), Zq_dot; Mv_dot, (Iy + Mq_dot)];

% The value of X_dot
f = [M_cal \ (A_cal + B_cal + C_cal + D_cal + E_cal + F_cal);
    -u0*sin(theta)+w*cos(theta); z; q];

% Computing Equilibrium Point
eqn = f([1:2]) == 0;
eqn = subs(eqn, [w, q], [0, 0]);

%Solving the equation for the equilibrium point
%S = solve(eqn, [w, q, theta, BAR])
S = solve(eqn, [theta, BAR])

% Equilibrium Points
q_e = 0
theta_e = double(S.theta)
w_e = u0 * tan(theta_e)
BAR_e = double(S.BAR)

% Computing the Jacobians
A = Jacobian(f, [w; q; z; int_z; theta])
B = Jacobian(f, BAR)

digits(4)
%Au = vpa(subs(A, [u0; theta; q; w], [u_val; theta_e; q_e; w_e]));
%Bu = vpa(subs(B, [u0; theta; q; w], [u_val; theta_e; q_e; w_e]));

Au = subs(A, [u0; theta; q; w], [u_val; theta_e; q_e; w_e]);
Bu = subs(B, [u0; theta; q; w], [u_val; theta_e; q_e; w_e]);
%disp(Au); disp(Bu)