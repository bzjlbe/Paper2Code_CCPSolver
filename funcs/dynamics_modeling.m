function [M_q, Q_q, pg_pq, v_n, v_t] = dynamics_modeling(q, dq, system_type, varargin)
% Dynamics modeling function
% Input parameters:
%   q: generalized coordinates
%   dq: generalized velocities
%   system_type: system type ('slider', 'Drod', or 'SP')
%   varargin: system-specific parameters
% Output parameters:
%   M_q: mass matrix
%   Q_q: generalized forces
%   pg_pq: Jacobian matrix for end point positions
%   v_n: normal velocities at contact points
%   v_t: tangential velocities at contact points

switch system_type
    case 'slider'
        % Slider system dynamics
        m = 1;
        g = 9.81;
        th_slope = varargin{1};
        
        % Mass matrix and generalized forces
        M_q = diag([m,m]);
        Q_q = -m*g*[sin(th_slope);cos(th_slope)];
        pg_pq = [1,0,0;0,0,1].';
        
        % Calculate end velocities
        v_end = pg_pq*dq;
        v_n = v_end(2);
        v_t = abs(v_end(1));
        
    case 'Drod'
        % Double rod system dynamics
        th = q(1);
        L = varargin{1};
        r = varargin{2};
        m = 0.402;
        J = m*L*L/12;
        g = 9.81;

        % Mass matrix and generalized forces
        M_q = diag([J,m,m]);
        Q_q = [0;0;-m*g];
        pg_pq = [-L*sin(th)/2-r*cos(th),0,-L*cos(th)/2+r*sin(th),-L*cos(th)/2-r*sin(th),0,L*sin(th)/2-r*cos(th);...
                 1,0,0,0,0,1;0,0,1,-1,0,0].';

        % Calculate end velocities
        v_end = pg_pq*dq;
        v_n = [v_end(3);v_end(6)];
        v_t = [(v_end(1)^2+v_end(2)^2)^0.5;(v_end(4)^2+v_end(5)^2)^0.5];
        
    case 'SP'
        % Slider-pendulum system dynamics
        th = q(3);
        dth = dq(3);
        m1 = varargin{1};
        m2 = varargin{2};
        nU = varargin{3};
        L = 2;
        J = m2*L^2/12;
        g = 10;

        % Mass matrix and generalized forces
        M_q = [m1+m2,0,m2*L/2*cos(th);0,m1+m2,m2*L/2*sin(th);m2*L/2*cos(th),0,J+m2*(L/2)^2];
        Q_q = [m2*L/2*dth^2*sin(th);-m2*L/2*dth^2*cos(th)-(m1+m2)*g;-m2*g*L/2*sin(th)];
        pg_pq = [1,0,0;0,0,1;0,0,0].';

        % Calculate end velocities
        v_end = reshape(pg_pq*dq, 3, nU).';
        v_t = (v_end(:,1).^2 + v_end(:,2).^2).^0.5;
        v_n = v_end(:,3);
        
    otherwise
        error('Unknown system type: %s', system_type);
end
end