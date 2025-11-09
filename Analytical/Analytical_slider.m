% --------------------------------------------------------------------- %
% This code is provided as supplementary information for the research article.
% The information about the article is as follows:
% article title: A prediction-correction method for spatial frictional contact problems in nonsmooth multibody dynamics
% --------------------------------------------------------------------- %

% Analytical Solution of Slider Motion
function output = Analytical_slider(miu,h,num,v0,th)
% miu: Friction coefficient
% h: The time step
% num: Number of steps
% v0: Slider tangential initial velocity
% th: Slope angle of the slope

a_G = -9.81;                            % Gravitational acceleration
a = a_G*sin(th) + miu*a_G*cos(th);      % Gravity
if a*v0 > 0
    a = -a;
end
t_mid = abs(v0/a);          % Slider stop moment
t_end = h*num;              % Simulation duration
hmi = h/1;
t1 = 0:hmi:t_mid;           % Sliding time discretization
t2 = t1(end)+h:hmi:t_end;   % Stationary time discretization

% Displacement
p1 = v0*t1 + 0.5*a*t1.^2;
p2 = -v0.^2/(2*a).*ones(size(t2));
% Velocity
v1 = v0+a*t1;
v2 = zeros(size(t2));
% Acceleration
d2q_real = [a(1)*ones(size(t1)),zeros(size(t2))];
% Contact forces
F_r = [[miu*a_G*cos(th);0;-a_G*cos(th)].*ones(size(t1)),[-a_G*sin(th);0;-a_G*cos(th)].*ones(size(t2))];

output.time = [t1,t2];
output.displacement = [p1,p2];
output.velocity = [v1,v2];
output.acceleration = d2q_real;
output.F_Analytical = F_r;