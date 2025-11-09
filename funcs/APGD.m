% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

function F_r = APGD(N,ss_r,tol,miu)
gd = 1e-6;
nC = size(N,1)/3;                   % Number of contact points
iter = 1;                           % Number of iterations
iter_max = 1000;                    % Maximum number of iterations
F_r = zeros(3*nC,1);                % Force/Impulse Iteration Initial Value
F_rback = ones(3*nC,1);             % fallback value
it_y = F_r;
it_th = 1;
it_L = norm(N*(F_r-F_rback))/norm(F_r-F_rback);
it_t = 1/it_L;                      % Iteration step
it_err_min = 1;                     % err
del_pr = [];

while iter < iter_max
    F_rlast = F_r;
    it_g = N*it_y + ss_r;
    % Projection del update force contact force F_r
    del = it_y - it_t*it_g;
    for i = 1:nC
        del_t = (del(3*(i-1)+1)^2 + del(3*(i-1)+2)^2)^0.5;
        if del_t < miu*del(3*i)                             % stick
            F_r(3*(i-1)+1:3*i,1) = del(3*(i-1)+1:3*i,1);
        elseif del_t < -del(3*i)/miu                        % separate
            F_r(3*(i-1)+1:3*i,1) = [0;0;0];
        else                                                % Slip
            F_r(3*(i-1)+1:3*i,1) = [miu*del(3*(i-1)+1)/del_t;miu*del(3*(i-1)+2)/del_t;1]*(miu*del_t+del(3*i))/(miu^2+1);
        end
    end
    % Check and adjust the step size it_t
    while 0.5*F_r.'*N*F_r+F_r.'*ss_r >= 0.5*it_y.'*N*it_y+it_y.'*ss_r + it_g.'*(F_r-it_y) + 0.5*it_L*norm(F_r-it_y)^2
        it_L = 2*it_L;
        it_t = 1/it_L;
        % Projection del update force contact force F_r
        del = it_y - it_t*it_g;
        for i = 1:nC
            del_t = (del(3*(i-1)+1)^2 + del(3*(i-1)+2)^2)^0.5;
            if del_t < miu*del(3*i)
                F_r(3*(i-1)+1:3*i,1) = del(3*(i-1)+1:3*i,1);
            elseif del_t < -del(3*i)/miu
                F_r(3*(i-1)+1:3*i,1) = [0;0;0];
            else
                F_r(3*(i-1)+1:3*i,1) = [miu*del(3*(i-1)+1)/del_t;miu*del(3*(i-1)+2)/del_t;1]*(miu*del_t+del(3*i))/(miu^2+1);
            end
        end
    end
    % Accelerated solution
    it_thlast = it_th;
    it_th = (-it_th^2+it_th*(it_th^2+4)^0.5)/2;
    it_ze = (1-it_thlast)*it_thlast/(it_thlast^2+it_th);
    it_y = F_r+it_ze*(F_r-F_rlast);
    % Calculate residuals
    del = F_r - gd*(N*F_r+ss_r);
    for i = 1:nC
        del_t = (del(3*(i-1)+1)^2 + del(3*(i-1)+2)^2)^0.5;
        if del_t < miu*del(3*i)
            del_pr(3*(i-1)+1:3*i,1) = del(3*(i-1)+1:3*i,1);
        elseif del_t < -del(3*i)/miu
            del_pr(3*(i-1)+1:3*i,1) = [0;0;0];
        else
            del_pr(3*(i-1)+1:3*i,1) = [miu*del(3*(i-1)+1)/del_t;miu*del(3*(i-1)+2)/del_t;1]*(miu*del_t+del(3*i))/(miu^2+1);
        end
    end
    it_err = norm((F_r-del_pr)/gd);
    % Fallback strategy to ensure residual monotonicity
    if it_err < it_err_min
        it_err_min = it_err;
        F_rback = F_r;
    end
    if it_err < tol
        break
    end
    iter = iter + 1;
    if it_g.'*(F_r-F_rlast) > 0
        it_y = F_r;
        it_th = 1;
    end
    it_L = 0.9*it_L;
    it_t = 1/it_L;
end
F_r = F_rback;