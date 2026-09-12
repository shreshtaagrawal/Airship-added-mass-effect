clear; clc;

% environment
p.rho   = 1.225;      % air density            [kg/m^3]
p.g     = 9.81;       % gravity                [m/s^2]

% geometry and mass
p.V_hull = 15.0;      % hull volume            [m^3]
p.L      = 8.0;       % hull length            [m]
p.m      = 18.0;      % vehicle mass           [kg]
p.z_b    = 0.30;      % CB height above CG     [m]
p.I_y    = 60.0;      % dry pitch inertia      [kg m^2]

% fluid inertia
p.k3      = 0.9;                        % transverse added-mass coefficient
p.m_added = p.k3 * p.rho * p.V_hull;    % added mass in heave     [kg]
p.I_added = 40.0;                       % added pitch inertia     [kg m^2]
p.m_eff   = p.m   + p.m_added;
p.I_eff   = p.I_y + p.I_added;

% derived 
p.B_force = p.rho * p.g * p.V_hull;     % buoyancy                [N]
p.k_pend  = p.B_force * p.z_b;          % pendulum stiffness      [N m/rad]
p.S_ref   = p.V_hull^(2/3);             % reference area          [m^2]

% aerodynamic coefficients
p.C_Zalpha =  2.00;   % normal force slope              [1/rad]
p.C_malpha =  0.02;   % Munk moment, POSITIVE = destabilising
p.C_mq     = -0.05;   % fin pitch damping

% actuator 
p.T   = 20.0;         % thrust magnitude                [N]
p.l_t = 1.0;          % thruster arm forward of CG      [m]



% Main block

[A_with, Bc, Bg] = build_model(3.0, p);

lambda = eig(A_with);
disp('Eigenvalues at U0 = 3 m/s:');
disp(lambda);

% independent pendulum prediction
wn_pend = sqrt(p.k_pend / p.I_eff);
fprintf('Pure pendulum frequency: %.4f rad/s\n', wn_pend);

% frequency of the oscillatory pair
osc = lambda(abs(imag(lambda)) > 1e-9);
if ~isempty(osc)
    fprintf('Model oscillatory |lambda|: %.4f rad/s\n', abs(osc(1)));
    fprintf('Damping ratio zeta      : %.4f\n', -real(osc(1))/abs(osc(1)));
end

%% ===== ADDED MASS SENSITIVITY STUDY =====

% Case 1: WITH added mass (physically correct)
[A_wet, ~, ~] = build_model(3.0, p);

% Case 2: WITHOUT added mass (dry)
p_dry        = p;           
p_dry.m_eff  = p.m;          % dry mass only
p_dry.I_eff  = p.I_y;        % dry inertia only
[A_dry, ~, ~] = build_model(3.0, p_dry);

% Extract the four metrics for each case
[wn_wet, T_wet, z_wet, tau_wet] = modes(A_wet);
[wn_dry, T_dry, z_dry, tau_dry] = modes(A_dry);

fprintf('\n%-24s %10s %10s %10s\n', 'Metric', 'WITH', 'WITHOUT', 'ratio');
fprintf('%-24s %10.4f %10.4f %10.3f\n', 'omega_n [rad/s]', wn_wet, wn_dry, wn_dry/wn_wet);
fprintf('%-24s %10.4f %10.4f %10.3f\n', 'period [s]',      T_wet,  T_dry,  T_dry/T_wet);
fprintf('%-24s %10.4f %10.4f %10.3f\n', 'damping zeta',    z_wet,  z_dry,  z_dry/z_wet);
fprintf('%-24s %10.4f %10.4f %10.3f\n', 'heave tau [s]',   tau_wet,tau_dry,tau_dry/tau_wet);




% Functions

function [A, Bc, Bg] = build_model(U0, p)

    qbar = 0.5 * p.rho * U0^2;

    % heave force derivatives
    Z_w     = -0.5 * p.rho * U0 * p.S_ref * p.C_Zalpha;
    Z_theta =  qbar * p.S_ref * p.C_Zalpha;
    Z_delta =  p.T;

    % pitch moment derivatives 
    M_w     = -0.5 * p.rho * U0 * p.S_ref * p.L * p.C_malpha;
    M_q     =  0.5 * p.rho * U0 * p.S_ref * p.L^2 * p.C_mq;
    M_munk  =  qbar * p.S_ref * p.L * p.C_malpha;
    M_theta = -p.k_pend + M_munk;
    M_delta =  p.T * p.l_t;

    % state matrix 
    A = [0,  1,               0,                   0;
         0,  Z_w/p.m_eff,     Z_theta/p.m_eff,     0;
         0,  0,               0,                   1;
         0,  M_w/p.I_eff,     M_theta/p.I_eff,     M_q/p.I_eff];

    Bc = [ 0; Z_delta/p.m_eff; 0; M_delta/p.I_eff];

    Bg = [ 0; -Z_w/p.m_eff; 0; -M_w/p.I_eff];

end

function [wn, T, zeta, tau] = modes(A)
    lam = eig(A);
    osc = lam(abs(imag(lam)) > 1e-9);     % the complex pair
    wn   = abs(osc(1));
    T    = 2*pi/abs(imag(osc(1)));
    zeta = -real(osc(1))/wn;
    
    real_poles = lam(abs(imag(lam)) < 1e-9 & abs(real(lam)) > 1e-9);
    tau  = -1/real(real_poles(1));        % heave time constant
end
