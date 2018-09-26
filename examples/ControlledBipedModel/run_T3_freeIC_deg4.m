%        theta = positive
% o     o     o
%  \    |--->/
%   \   ^   / l
%    \  F  /
%     \ | /
%      \|/
% ============

polysin = @(ang) ang - ang.^3/6;
polycos = @(ang) 1 - ang.^2/2;

%-------------------------------------------------------------------------%
%------------------------ All Physical Parameters ------------------------%
%-------------------------------------------------------------------------%
params = struct();

params.m        = 1;            % mass
params.g        = 0.2;          % gravitational acceleration
params.l0       = 0.5;          % leg length at touch-down
params.lmax     = 1;            % maximum leg length

params.alpha    =-pi/5;         % leg angle in flight phase = -30 degrees
params.umax     = 1;            % upper bound of input

%-------------------------------------------------------------------------%
%-------------- Domains (Defined by Upper and Lower Bounds) --------------%
%-------------------------------------------------------------------------%
al = params.alpha;
yR_lo = params.l0 * polycos(al);        % touch-down height
yR_hi = yR_lo + 0.04;
params.yR_lo = yR_lo;
params.yR_hi = yR_hi;
lmax = params.lmax;


% Mode 1: stance, y <= yR
% state = ( l, l_dot, theta, theta_dot )
params.domain{1} =...
        [ 0.1, lmax;        % l         - leg length
         -0.6, 0.6;         % l_dot     - time derivative of l
         -0.7, 1.2;         % theta     - leg angle
         -0.1, 2 ];         % theta_dot - time derivative of theta

% Mode 2: stance, y >= yR
% state = ( l, l_dot, theta, theta_dot )
params.domain{2} =...
        [ 0.1, lmax;        % l         - leg length
         -0.6, 0.6;         % l_dot     - time derivative of l
         -0.7, 1.2;         % theta     - leg angle
         -0.1, 2 ];         % theta_dot - time derivative of theta

%-------------------------------------------------------------------------%
%-------------------------- Parameters for OCP ---------------------------%
%-------------------------------------------------------------------------%
d = 4;              % degree of relaxation
T = 3;              % time horizon
nmodes = 2;         % number of modes

% Solver options
options.freeFinalTime = 0;      % fixed terminal time
options.withInputs = 1;         % control extraction?
options.svd_eps = 1e4;          % svd threshould for moment matrices

%-------------------------------------------------------------------------%
%---------------------------- Construct OCP ------------------------------%
%-------------------------------------------------------------------------%
% Define variables
t = msspoly( 't', 1 );
x = cell( nmodes, 1 );
u = cell( nmodes, 1 );
f = cell( nmodes, 1 );
g = cell( nmodes, 1 );
x0 = cell( nmodes, 1 );
hX0 = cell( nmodes, 1 );
hX = cell( nmodes, 1 );
hU = cell( nmodes, 1 );
hXT = cell( nmodes, 1 );
sX = cell( nmodes, nmodes );
R = cell( nmodes, nmodes );
h = cell( nmodes, 1 );
H = cell( nmodes, 1 );
c = cell( nmodes, 1 );

xvar = msspoly('x', 4);
uvar = msspoly('u', 1);

x{1} = xvar;
u{1} = uvar;
x{2} = xvar;
u{2} = uvar;

% Dynamics
for i = 1 : 2
    f{i} = T * Swing_f_poly( x{i}, params );
    g{i} = T * Swing_g_poly( x{i}, params );
end

% Suppports, Reset Maps, and Cost Functions
domain  = params.domain;
l0      = params.l0;
umax    = params.umax;
al      = params.alpha;

y = xvar(1) * polycos(xvar(3));             % y = l * cos(theta)

d_des = 1;        % desired step length

% Mode 1 : Stance (y <= yR_hi)
hX{1} = [ domain{1}(:,2) - xvar;            % domain
          xvar - domain{1}(:,1);
          yR_hi - y ];
hU{1} = uvar * (umax - uvar);
R{1,2} = xvar;                              % reset map
sX{1,2} = ...                               % guard
        [ y - yR_hi;             	% y = yR
          yR_hi - y;
          hX{1} ];                  % G \subset X
c{1,2} = msspoly(0);
h{1} = uvar^2 * T;
H{1} = msspoly(0);

% Mode 2 : Stance (y >= yR_lo)
hX{2} = [ domain{1}(:,2) - xvar;            % domain
          xvar - domain{1}(:,1);
          y - yR_lo ];
hU{2} = uvar * (umax - uvar);
R{2,1} = Reset_poly(xvar,params);           % reset map
sX{2,1} = ...                               % guard
        [ y - yR_lo;                % y = yR
          yR_lo - y;
          hX{2} ];              	% G \subset X
c{2,1} = 10 * (xvar(1) * polysin(xvar(3)) + l0 * polysin(-al) - d_des) ^ 2;        % step length = l * sin(theta) + l0 * sin(-alpha)
h{2} = uvar^2 * T;
H{2} = msspoly(0);


% Initial condition and Target Set
x0_ub = [0.35; 0; 0; 1.2];
x0_lb = [0.35; 0; 0; 0.7];
hX0{1} = [ x0_ub - xvar
           xvar - x0_lb ];

% Target set is the entire X{1}
hXT{1} = hX{1};

% Solve
[out] = HybridOCPDualSolver_switching_IC(t,x,u,f,g,hX,hU,sX,R,hX0,hXT,h,H,c,d,options);

disp(['p* = ', num2str( out.pval )]);

fprintf('Finished solving freeIC, T%d, deg %d, F approx 3.\n', T, d);
save(['Data/Result_T',num2str(T),'_freeIC_deg',num2str(d)]);