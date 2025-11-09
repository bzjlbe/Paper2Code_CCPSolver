function output = PM_nonsm_slider(miu, rho, h, num, v0, th_slope, sw_on)
% Main file for slider system
% Input parameters:
%   miu: friction coefficient
%   rho: generalized_alpha method parameter
%   h: time step size
%   num: number of steps
%   v0: initial tangential velocity of slider
%   th_slope: slope angle
%   sw_on: algorithm switch (0 for CCP_Wang, 1 for proposed method)

tic
% Basic parameters
nq = 2; nU = 1;

% Initial conditions
total_steps = num + 1;
[q, dq, dq_sm, dq_imp, d2q, zeta_0, zeta] = deal(zeros(nq, total_steps));
[F_end, P_end] = deal(zeros(3*nU, total_steps));
[v_n, v_t, dmvt] = deal(zeros(nU, 1));
dq(:,1) = [v0;0];
dq_sm(:,1) = [v0;0];
Fyi = zeros(nU,num+1);


% Configuration parameters
params = config(rho, h, num, 0.8, 1e-8);

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
    [nC, nI, J_Fc_F, J_Fi_F, J_Fa_F, Fyi] = contact_detection(q(:,n+1), Fyi, n, nU, 'slider');

    % Dynamics solving
    [q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t, dmvt_n] = dynamics_solver(...
        q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t, params, n, miu, J_Fc_F, J_Fi_F,...
        J_Fa_F, nC, nI, params.e_imp, params.tol, sw_on, 'slider', th_slope);
    dmvt(:,n+1) = dmvt_n;
end

toc

% Output results
output.miu = miu;
output.time = params.time;
output.q = q;
output.dq = dq;
output.d2q = d2q;
output.th = q(1,:);
output.x = q(2,:);
output.dth = dq(1,:);
output.dx = dq(2,:);
output.zeta_0 = zeta_0;
output.F_end = F_end;
output.dmvt = dmvt;
end