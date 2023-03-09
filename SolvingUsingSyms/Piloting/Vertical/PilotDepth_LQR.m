function out = PilotDepth_LQR(in, memory, parameters, veh_model)
% Depth piloting
%coder.extrinsic('lqrd');
coder.extrinsic('computeABCMatrices');

%Q = diag([1e-6, 1,100,10,0.1]);
%R = 1e9;

Q = 1 * diag([10,1,100,10,100]);
R = 1e9;

%% Loading Parameters

u0 = in.u_ms;

w_e = 0; q_e = 0; theta_e = 0; BAR0 = 0; K = zeros(1,5); 
[w_e, q_e, theta_e, BAR0, K] = computeABCMatrices(in, veh_model, u0, Q, R);

%%
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = w_e - in.w_ms;
delta_q = q_e - in.q_rads;
delta_theta = DiffAngle(theta_e,in.theta_rad);

BARc_rad = BAR0 - delta_w.*K(1) - delta_q.*K(2) - delta_z.*K(3)...
   - delta_theta.*K(4)  - memory.int_z.*K(5) ;

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
