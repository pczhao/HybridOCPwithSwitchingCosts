% Solve biped model HOCP with GPOPS-II

clear;
clc;

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
params.umax     = 1;          % upper bound of input

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
%--------------- Provide All Physical Data for Problem -------------------%
%-------------------------------------------------------------------------%
T = 6;
d_des = 1;
nphases = 5;
x0 = [ 0.35, 0, 0, 0.85 ];
x0_ub = [0.35; 0; 0; 1.2]';
x0_lb = [0.35; 0; 0; 0.7]';

auxdata = struct;
auxdata.params = params;
auxdata.l0 = params.l0;
auxdata.yR_lo = params.yR_lo;
auxdata.yR_hi = params.yR_hi;
auxdata.nphases = nphases; 
auxdata.T = T;
auxdata.d_des = d_des;
auxdata.alpha = params.alpha;

%-------------------------------------------------------------------------%
%----------------- Generate Initial Guess for Problem --------------------%
%-------------------------------------------------------------------------%
clear guess bounds mesh
dt = (T/1.5) / nphases;
for iphase = 1 : nphases
    current_mode = mod( iphase+1, 2 ) + 1;
    domain = params.domain{current_mode};
    
    guess.phase(iphase).time = [ dt*(iphase-1); dt*iphase ];
    guess.phase(iphase).state = zeros(2,4);
    guess.phase(iphase).control = [ 0; 0 ];
    guess.phase(iphase).integral = 0;
end

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%

t0 = 0;

% Mode:  2 -> 1  (In spotless def)
% Phase: 1 -> 2  (In gpops def)

for iphase = 1 : nphases
    idx = mod( iphase+1, 2 ) + 1;
    domain = params.domain{idx};
    bounds.phase(iphase).initialtime.lower = t0;
    bounds.phase(iphase).finaltime.upper = T;
    if iphase == 1
        bounds.phase(iphase).initialtime.upper = t0;
    else
        bounds.phase(iphase).initialtime.upper = T; 
    end
    if iphase == nphases
        bounds.phase(iphase).finaltime.lower = T;
    else
        bounds.phase(iphase).finaltime.lower = t0;
    end
    if iphase == 1
        bounds.phase(iphase).initialstate.lower = x0_lb;
        bounds.phase(iphase).initialstate.upper = x0_ub;
    else
        bounds.phase(iphase).initialstate.lower = domain(:,1)';
        bounds.phase(iphase).initialstate.upper = domain(:,2)';
    end
    bounds.phase(iphase).state.lower = domain(:,1)';
    bounds.phase(iphase).state.upper = domain(:,2)';
    bounds.phase(iphase).finalstate.lower = domain(:,1)';
    bounds.phase(iphase).finalstate.upper = domain(:,2)';
    bounds.phase(iphase).control.lower = 0;
    bounds.phase(iphase).control.upper = 1;
    bounds.phase(iphase).integral.lower = -100000;
    bounds.phase(iphase).integral.upper = 100000;
    
    if (idx == 1)
        % y \in [0, yR_hi]
        bounds.phase(iphase).path.lower = 0;
        bounds.phase(iphase).path.upper = params.yR_hi;
    else
        % y \in [yR_lo, lmax]
        bounds.phase(iphase).path.lower = params.yR_lo;
        bounds.phase(iphase).path.upper = params.lmax;
    end
    bounds.phase(iphase).duration.lower = 0.05;
    bounds.phase(iphase).duration.upper = T;
end

for iphase = 1 : nphases-1
    idx = mod( iphase+1, 2 ) + 1;
    switch idx
        case 1              % {y<=yR_hi} -> {y>=yR_lo}
            bounds.eventgroup(iphase).lower = zeros(1,6);
            bounds.eventgroup(iphase).upper = zeros(1,6);
            
        case 2              % {y>=yR_lo} -> {y<=yR_lo}
            bounds.eventgroup(iphase).lower = zeros(1,6);
            bounds.eventgroup(iphase).upper = zeros(1,6);
            
    end
end

% Terminal condition
iphase = nphases;
bounds.eventgroup(iphase).lower = 0;
bounds.eventgroup(iphase).upper = 1000;

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
% mesh.maxiterations   = 45;
mesh.tolerance       = 1e-7;
for i = 1 : nphases
    mesh.phase(i).colpoints = 10 * ones(1,100);
    mesh.phase(i).fraction = 0.01 * ones(1,100);
end


%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'BipedModel';
setup.functions.continuous           = @BipedModelContinuous;
setup.functions.endpoint             = @BipedModelEndpoint;
setup.displaylevel                   = 2;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.auxdata                        = auxdata;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.maxiterations = 2000;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-7;
% setup.derivatives.supplier           = 'adigator';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
% setup.scales.method                  = 'automatic-hybridUpdate';


%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc

disp(['Gpops p* = ', num2str(output.result.objective)]);

%-------------------------------------------------------------------------%
%-------------------------------- PLot -----------------------------------%
%-------------------------------------------------------------------------%
%%
state_hist_gpops = [];
t_hist_gpops = [];
control_hist_gpops = [];
xoffset = [];
phase_hist_gpops = [];
previous_x = 0;
for iphase = 1 : nphases
    t_hist_gpops = [ t_hist_gpops; output.result.solution.phase(iphase).time ];
    state_hist_gpops = [ state_hist_gpops; output.result.solution.phase(iphase).state ];
    control_hist_gpops = [ control_hist_gpops; output.result.solution.phase(iphase).control ];
    xoffset = [ xoffset; previous_x * ones(length(output.result.solution.phase(iphase).time),1) ];
    if (mod(iphase,2) == 0)
        xf1 = state_hist_gpops( end, : );
        previous_x = previous_x + xf1(1) * polysin(xf1(3)) + params.l0 * polysin(-params.alpha);
    end
    phase_hist_gpops = [ phase_hist_gpops; iphase * ones(length( output.result.solution.phase(iphase).time ), 1) ];
end
l_hist = state_hist_gpops( :, 1 );
ldot_hist = state_hist_gpops( :, 2 );
theta_hist = state_hist_gpops( :, 3 );
thetadot_hist = state_hist_gpops( :, 4 );
x_hist = l_hist .* polysin( theta_hist ) + xoffset;
y_hist = l_hist .* polycos( theta_hist );

save(['Result_T',num2str(T),'_freeIC_gpops']);


figure;
hold on;
plot(x_hist, y_hist);
title('X-Y');
