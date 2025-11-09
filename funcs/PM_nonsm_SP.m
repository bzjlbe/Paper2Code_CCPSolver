function output = PM_nonsm_SP(miu, th0, dx0, m1, m2, rho, h, num, sw_on)
% Main file for slider-pendulum system
% Input parameters:
%   miu: friction coefficient
%   th0: initial rod angle
%   dx0: initial horizontal velocity of slider
%   m1: slider mass
%   m2: rod mass
%   rho: generalized_alpha method parameter
%   h: time step size
%   num: number of steps
%   sw_on: algorithm switch (0 for CCP_Wang, 1 for proposed method)

tic
% Basic parameters
nq = 3; nU = 1;

% Initial conditions
x0 = 0; y0 = 0; dy0 = 0; dth0 = 0;
total_steps = num + 1;
[q, dq, dq_sm, dq_imp, d2q, zeta_0, zeta] = deal(zeros(nq, total_steps));
[F_end, P_end] = deal(zeros(3*nU, total_steps));
[v_n, v_t] = deal(zeros(nU, 1));
Fyi = zeros(nU, total_steps);
q(:,1) = [x0;y0;th0];dq(:,1) = [dx0;dy0;dth0];dq_sm(:,1) = [dx0;dy0;dth0];

% Configuration parameters
params = config(rho, h, num, 0.8, 1e-6);

for n = 0:num
    % Generalized-alpha algorithm update
    if n ~= 0
        zeta_0(:,n+1) = zeta(:,n) + (1-params.alpha_f)/(1-params.alpha_m)*d2q(:,n);
        zeta(:,n+1) = (params.alpha_f*d2q(:,n) - params.alpha_m*zeta_0(:,n+1))/(1-params.alpha_m);
        q(:,n+1) = q(:,n) + h*dq(:,n) + h*h*(1-2*params.eta)*zeta_0(:,n+1)/2 + h*h*params.eta*zeta(:,n+1);
        dq_sm(:,n+1) = dq(:,n) + h*(1-params.gama)*zeta_0(:,n+1) + h*params.gama*zeta(:,n+1);
        dq(:,n+1) = dq_sm(:,n+1);
    end
    
    % Contact detection
    [nC, nI, J_Fc_F, J_Fi_F, J_Fa_F, Fyi] = contact_detection(q(:,n+1), Fyi, n, nU, 'SP');

    % Dynamics solving
    [q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t] = dynamics_solver(...
        q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t, params, n, miu, J_Fc_F, J_Fi_F,...
        J_Fa_F, nC, nI, params.e_imp, params.tol, sw_on, 'SP', m1, m2, nU);
end

toc

% Output results
output.time = params.time;
output.th = q(3,:);
output.x = q(1,:);
output.y = q(2,:);
output.dth = dq(3,:);
output.dx = dq(1,:);
output.dy = dq(2,:);
output.d2th = d2q(3,:);
output.d2x = d2q(1,:);
output.d2y = d2q(2,:);
output.Fx = F_end(1,:);
output.Fz = F_end(3,:);
output.Px = P_end(1,:);
output.Pz = P_end(3,:);
output.zeta_0 = zeta_0;
end