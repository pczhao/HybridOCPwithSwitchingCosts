% draw biped model

% Generate sample trajectory

clear;
close all;
l0 = 0.5;           % normal
la = 0.04;           % actuator length
w = 0.04;
h = 0.07;

theta = pi/7;

h_fig = figure;
hold on;
%% Plot left
Opos = [ 0, 0 ];

Mpos = Opos + l0 * [ sin(theta), cos(theta) ];
Fpos = Opos + la * [ -sin(theta), cos(theta) ];

RO = [  cos(theta), sin(theta);
       -sin(theta), cos(theta) ];


% Below actuator
lb = (l0 - la) / 2;
lineB = [ 0, 0
          0, lb ]';
tmp = RO * lineB;
plot(tmp(1,:), tmp(2,:), 'k', 'LineWidth', 2);

axis equal

% Above actuator
lineA = [ 0, l0
          0, l0-lb  ]';
tmp = RO * lineA;
plot(tmp(1,:), tmp(2,:), 'k', 'LineWidth', 2);

% Actuator
line1 = [ -w/2, l0-lb
           w/2, l0-lb ]';
line2 = [ -w/2, lb+h
          -w/2, lb
           w/2, lb
           w/2, lb+h ]';
tmp = RO * line1;
plot(tmp(1,:), tmp(2,:), 'k', 'LineWidth', 2);
tmp = RO * line2;
plot(tmp(1,:), tmp(2,:), 'k', 'LineWidth', 2);

%% Plot right
al = -pi/5;

RM = [  cos(al), sin(al), Mpos(1);
       -sin(al), cos(al), Mpos(2);
              0,       0,       1 ];

lineS = [ 0, 0, 1
          0, -l0, 1]';
tmp = RM * lineS;
plot(tmp(1,:), tmp(2,:), 'Color', [0.6,0.6,0.6], 'LineWidth', 2);

% Mass
th = 0:pi/50:2*pi;
xvec = Mpos(1) + 0.03*cos(th);
yvec = Mpos(2) + 0.03*sin(th);
patch(xvec, yvec, 1, 'FaceColor', [0.8,0.8,0.8], 'LineWidth',2);


%% Draw ground
% Draw the ground. It reaches from -2.5 to +6.5.
            h   = 0.005; % Height of the bar at the top
            n   = 5000;  % Number of diagonal stripes in the shaded area
            s   = 0.05; % Spacing of the stripes
            w   = 0.005; % Width of the stripes
            ext = 0.02;  % Length of the stripes
    
            % Create vertices by shifting a predefined pattern 'n' times to the right:
            v = [     -50,0;
                repmat([     0,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
                repmat([  -ext,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
                repmat([-ext+w,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
                repmat([     w,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
                -50+s*n+w,0];
            % Connect to faces:
            f = [1,2,4*n+1,4*n+2;
                 repmat([0,n,2*n,3*n],n,1) + repmat((1:n)',1,4)+1];

            
            % Define some arbitrary states:
            x        = 1;
            y        = 1.2;
            l_leg    = 1;
            phi_body = 0;
                        
            
            % Draw Ground
            vert_x_out = [-15 100 100 -15];
            vert_y_out = [0 0 -20 -20];
%             patch(vert_x_out, vert_y_out,'white');   
            patch('faces', f, 'vertices', v); 

%% Auxiliary
plot([0,0],[0,0.2], 'k-','LineWidth',1);
Spos = Mpos + l0 * [sin(-al), -cos(al)];
plot([Spos(1), Spos(1)], [Spos(2), Spos(2)+0.2], 'k--','LineWidth', 1);

plot([0.6, 0.6], [0, 0.02], 'k-', 'LineWidth', 1);
%% Size
axis equal
% xlim([0,0.6]);
xlim([-0.1,0.7]);
% xlim([0.3,0.6]);
ylim([-0.12,0.57]);
axis off

set(h_fig,'PaperSize',[4.5 3.5]); %set the paper size to what you want  
h_fig.PaperPositionMode = 'auto';
% fig_pos = h_fig.PaperPosition;
% h_fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(h_fig,'ModelDef','-dpdf') % then print it
