function FP_r = cone_complementarity_solver(N, ss_r1, miu, v_t, v_n, n, e_imp, tol, sw_on, J_Fc_F, J_Fi_F, J_Fa_F)
% Cone complementarity problem solver
% Input parameters:
%   N: system matrix
%   ss_r1: intermediate item
%   miu: friction coefficient
%   v_t, v_n: tangential and normal velocities
%   n: current time step
%   e_imp: restitution coefficient
%   tol: tolerance
%   sw_on: algorithm switch (0 for CCP_Wang, 1 for proposed method)
%   J_Fc_F, J_Fi_F, J_Fa_F: selection matrices

b_con = J_Fc_F*reshape([zeros(size(v_t,1),2), miu*v_t(:,n+1)].', 3*size(v_t,1), 1);
if n ~= 0
    b_imp = J_Fi_F*reshape([zeros(size(v_t,1),2), miu*v_t(:,n+1) + e_imp*min(v_n(:,n),0)].', 3*size(v_t,1), 1);
else
    b_imp = J_Fi_F*reshape([zeros(size(v_t,1),2), miu*v_t(:,n+1)].', 3*size(v_t,1), 1);
end
b_all = [b_con; b_imp];
ss_r = ss_r1 + b_all;
FP_r = APGD(N, ss_r, tol, miu);

% Middle iteration (executed only when sw_on=1)
for k_me = 1:100*sw_on
    b_last = b_all;
    v_end = reshape(J_Fa_F.'*(N*FP_r + ss_r1), 3, size(v_t,1)).';
    v_t_new = (v_end(:,1).^2 + v_end(:,2).^2).^0.5;

    % Update friction cone boundaries
    b_con = J_Fc_F*reshape([zeros(size(v_t_new,1),2), miu*v_t_new].', 3*size(v_t_new,1), 1);
    if n ~= 0
        b_imp = J_Fi_F*reshape([zeros(size(v_t_new,1),2), miu*v_t_new + e_imp*min(v_n(:,n),0)].', 3*size(v_t_new,1), 1);
    else
        b_imp = J_Fi_F*reshape([zeros(size(v_t_new,1),2), miu*v_t_new].', 3*size(v_t_new,1), 1);
    end
    b_all = [b_con; b_imp];
    
    if norm(b_all - b_last) < tol
        break
    end
    
    ss_r = ss_r1 + b_all;
    FP_r = APGD(N, ss_r, tol, miu);
end
end