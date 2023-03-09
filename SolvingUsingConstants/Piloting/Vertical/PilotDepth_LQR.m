function out = PilotDepth_LQR(in, memory, parameters)

%% Loading Parameters
w_e = 0.0138;
q_e = 0;
theta_e = 0.0069;
BAR_e = -0.0837;
%K = [-0.0051    0.0001   -0.0022   -0.0032   -0.0000];
%K = [-0.0447   -0.0927    0.0022   -0.1341    0.0000];
K = [-0.0887   -0.0170   -0.0024   -0.0060   -0.0000];

%%
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = w_e - in.w_ms;
delta_q = q_e - in.q_rads;
delta_theta = DiffAngle(theta_e,in.theta_rad);

BARc_rad = BAR_e - delta_w.*K(1) - delta_q.*K(2) - delta_z.*K(3)...
   - delta_theta.*K(4)  - memory.int_z.*K(5) ;

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
