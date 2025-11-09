function [nC, nI, J_Fc_F, J_Fi_F, J_Fa_F, Fyi] = contact_detection(q, Fyi, n, nU, system_type, varargin)
% Contact detection function
% Input parameters:
%   q: generalized coordinates
%   Fyi: constraint switch history
%   n: current time step
%   nU: number of potential contact points
%   system_type: system type ('slider', 'Drod', 'SP')
%   varargin: additional parameters (system dependent)

Fyi_current = zeros(nU,1);
UnitnU = eye(3*nU);

switch system_type
    case 'slider'
        % Slider system contact detection
        Fyi_current(:,1) = (q(2) <= 0e-6);
        
    case 'Drod'
        % Double rod system contact detection
        [L, r] = varargin{1:2};% rod length % rod cross-section radius
        [th, x, z] = deal(q(1), q(2), q(3));

        % Calculate end positions
        r_end(1:3,1) = [x;0;z] + [L*cos(th)/2-r*sin(th);0;-L*sin(th)/2-r*cos(th)];
        r_end(4:6,1) = [-z;0;x] + [-L*sin(th)/2+r*cos(th);0; -L*cos(th)/2-r*sin(th)];
        Fyi_current(:,1) = ([r_end(3,1);r_end(6,1)] <= 1e-4);
        
    case 'SP'
        % Slider-pendulum system contact detection
        Fyi_current(:,1) = (q(2) <= 1e-4);
end

% Update constraint history
Fyi(:,n+1) = Fyi_current;

% Determine contact types
if n == 0
    gr_C = 2 * Fyi_current;
    nC = sum(gr_C == 2);
    nI = 0;
else
    gr_C = 2 * Fyi_current + Fyi(:,n);
    nC = sum(gr_C == 3);    % number of continuous contact points
    nI = sum(gr_C == 2);    % number of non-smooth contact points
end

% nA = nC + nI;

% Construct selection matrices
Fyi_nA = find(gr_C >= 2);
Fyi_3nA = reshape([3*Fyi_nA-2,3*Fyi_nA-1,3*Fyi_nA].',3*numel(Fyi_nA),1);
J_Fa_F = UnitnU(Fyi_3nA,:);
if nI == 0
    J_Fc_F = J_Fa_F;                % contact point selection matrix
    J_Fi_F = double.empty(0,3*nU);  % non-smooth point selection matrix
else
    J_Fc_F = double.empty(0,3*nU);
    J_Fi_F = J_Fa_F;
end
end