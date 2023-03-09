function [w_e, q_e, theta_e, BAR_e, K] = computeABCMatrices(in, veh_model, u_val, Q, R)

if ~isfield(veh_model, 'matrices')

    %A = double(subs(A_mat, u0, u_val));
    %B = double(subs(B_mat, u0, u_val));

    %A = double(A_mat(u_val));
    %B = double(B_mat(u_val));

    rho = veh_model.Environment.Rho;
    % nu = veh_model.Environment.nu;
    Sref = veh_model.Mecanic.surface_reference;
    Lref = veh_model.Mecanic.length_reference;

    CMuw = veh_model.Hydrodynamic.CMuw;
    CMuq = veh_model.Hydrodynamic.CMuq;

    CZuw = veh_model.Hydrodynamic.CZuw;
    CZuq = veh_model.Hydrodynamic.CZuq;

    CZ0 = veh_model.Hydrodynamic.CZ0;
    CM0 = veh_model.Hydrodynamic.CM0;

    m = veh_model.Mecanic.mass;
    Zw_dot = veh_model.Mecanic.added_water_matrix(3, 3);
    Zq_dot = veh_model.Mecanic.added_water_matrix(3, 5);
    Mv_dot = veh_model.Mecanic.added_water_matrix(5, 2);
    Iy = veh_model.Mecanic.inertia_matrix(2, 2);
    Mq_dot = veh_model.Mecanic.added_water_matrix(5, 5);
    zg = veh_model.Mecanic.center_of_gravity_veh_ref(3);
    zb = veh_model.Mecanic.center_of_buoyancy_veh_ref(3);

    CZ = veh_model.XRearHelms.CZ;
    CM = veh_model.XRearHelms.CM;

    gravity = 9.81;
    V = veh_model.Mecanic.volume;
    
    W = m * gravity;
    B = V*rho*gravity;
    
    k = 0.5*rho*Sref;
    
    %% Etats Equilibres
    
    syms u0 w q z int_z theta BARx
    
    A_cal = u0 * [0; m*w];
    B_cal = [0, m*(zg*q + u0); -m*(zg*q + u0), 0] * [w; q];
    C_cal = k * [CZ0*abs(u0); Lref * CM0 * abs(u0)] * u0;
    D_cal = k * [CZuw * u0, CZuq*Lref*u0; CMuw*Lref*u0, Lref^2*CMuq*u0]*[w; q];
    E_cal = [(W - B)*cos(theta); -(zg*W - zb*B)*sin(theta)];
    F_cal = k * [CZ*BARx*u0* abs(u0); Lref*CM*BARx*u0*abs(u0)];
    M_cal = [(m+Zw_dot), Zq_dot; Mv_dot, (Iy + Mq_dot)];
    
    % The value of X_dot
    f = [M_cal \ (A_cal + B_cal + C_cal + D_cal + E_cal + F_cal);
        w; q; z];
    
    % Computing Equilibrium Point
    eqn = f(1:4) == 0;
    eqn = subs(eqn, u0, u_val);
    
    %Solving the equation for the equilibrium point
    S = vpasolve(eqn, [w, q, theta, BARx]);
    
    % Equilibrium Points
    w_e = double(S.w);
    q_e = double(S.q);
    theta_e = double(S.theta);
    BAR_e = double(S.BARx);
    
    %% A, B and C matrices
    
    % Computing the Jacobians
    A_mat = Jacobian(f, [w; q; z; theta; int_z]);
    B_mat = Jacobian(f, BARx);
    
    digits(4)
    
    A_mat = subs(A_mat, [theta; q; w; BARx], [theta_e; q_e; w_e; BAR_e]);
    B_mat = subs(B_mat, [theta; q; w; BARx], [theta_e; q_e; w_e; BAR_e]);

    C = eye(5);
    
    %veh_model.matrices.A = matlabFunction(A_mat);
    %veh_model.matrices.B = matlabFunction(B_mat);
    
    veh_model.matrices.C = C;

    veh_model.matrices.w_e = w_e;
    veh_model.matrices.q_e = q_e;
    veh_model.matrices.theta_e = theta_e;
    veh_model.matrices.BAR_e = BAR_e;

    A = double(subs(A_mat, u0, u_val));
    B = double(subs(B_mat, u0, u_val));

    veh_model.matrices.A = A;
    veh_model.matrices.B = B;

    Qy = C'*Q*C;

    K = zeros(1,5);
    K = lqrd(A,B,Qy,R,in.delta_time_s);

    veh_model.matrices.K = K;

else

    w_e = veh_model.matrices.w_e;
    q_e = veh_model.matrices.q_e;
    theta_e = veh_model.matrices.theta_e;
    BAR_e = veh_model.matrices.BAR_e;

    K = veh_model.matrices.K;

end

end

