function [q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t, dmvt] = dynamics_solver(q, dq, d2q, dq_imp, F_end, P_end, zeta_0, zeta, v_n, v_t, ...
    params, n, miu, J_Fc_F, J_Fi_F, J_Fa_F, nC, nI, e_imp, tol, sw_on, system_type, varargin)
% Dynamics solver for one time step with outer iteration loop
% Outer iteration loop
while(1)
    % Dynamics modeling - general call to dynamics modeling function
    [M_q, Q_q, pg_pq, vn, vt] = dynamics_modeling(q(:,n+1), dq(:,n+1), system_type, varargin{:});
    
    % Store calculated velocities
    v_n(:,n+1) = vn;
    v_t(:,n+1) = vt;

    % Calculate residuals
    F_sm_Uni = J_Fc_F*F_end(:,n+1);
    P_nsm_Uni = J_Fi_F*P_end(:,n+1);

    res_C = M_q*d2q(:,n+1) - (J_Fc_F*pg_pq).'*F_sm_Uni - Q_q;
    res_I = M_q*dq_imp(:,n+1) - (J_Fi_F*pg_pq).'*P_nsm_Uni;
    res = [res_C; res_I];

    if norm(res,2) < tol
        break
    end

    % Store previous values
    F_end0 = F_end(:,n+1);
    P_end0 = P_end(:,n+1);

    % Construct system matrices
    J_q = [params.gama_t*M_q, zeros(size(M_q)); -M_q, M_q];
    D_con = J_Fc_F*pg_pq;
    D_imp = J_Fi_F*pg_pq;
    D_all = [D_con; D_imp];
    B_all = blkdiag(D_con.', D_imp.');
    trC_all = [zeros(3*(nC+nI), size(M_q,2)), D_all];
    invJ_q = J_q^(-1);

    N = trC_all*invJ_q*B_all;

    % Construct friction cone boundaries
    ss_r1 = D_all*dq(:,n+1) - trC_all*invJ_q*(B_all*blkdiag(J_Fc_F,J_Fi_F)*[F_end0;P_end0] + res);

    mvt1 = miu*v_t(:,n+1);
    % Solve cone complementarity problem
    FP_r = cone_complementarity_solver(N, ss_r1, miu, v_t, v_n, n, e_imp, tol, sw_on, J_Fc_F, J_Fi_F, J_Fa_F);

    v_end = reshape(J_Fa_F.'*(N*FP_r + ss_r1), 3, size(v_t,1)).';
    v_t_new = (v_end(:,1).^2 + v_end(:,2).^2).^0.5;
    mvt2 = miu*v_t_new;
    dmvt = mvt1 - mvt2;

    % Update contact forces/impulses
    if nI == 0
        F_end(:,n+1) = J_Fa_F.'*FP_r;
    else
        P_end(:,n+1) = J_Fa_F.'*FP_r;
    end

    % Update state variables
    del_F_alpha = [J_Fc_F*(F_end(:,n+1) - F_end0); J_Fi_F*(P_end(:,n+1) - P_end0)];
    del_dx = invJ_q*(B_all*del_F_alpha - res);

    if n ~= 0
        del_dq_sm = del_dx(1:size(M_q,2),:);
        del_dq = del_dx(size(M_q,2)+1:end,:);
        del_dq_imp = del_dq - del_dq_sm;

        d2q(:,n+1) = d2q(:,n+1) + del_dq_sm*params.gama_t;
        dq(:,n+1) = dq(:,n+1) + del_dq;
        dq_imp(:,n+1) = dq_imp(:,n+1) + del_dq_imp;
        q(:,n+1) = q(:,n+1) + del_dq_sm*params.eta_t + (params.h/2)*del_dq_imp;
    else
        del_dq_sm = del_dx(1:size(M_q,2),:);
        d2q(:,n+1) = d2q(:,n+1) + del_dq_sm*params.gama_t;
        zeta_0(:,n+1) = d2q(:,n+1);
        zeta(:,n+1) = (params.alpha_f*d2q(:,n+1) - params.alpha_m*zeta_0(:,n+1))/(1-params.alpha_m);
    end
end
end