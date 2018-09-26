function plotdata_freeIC( h_axis, filename, xs, mycolor )
if (nargin < 2) || (isempty(filename))
    load('Result_Fpoly3_T7_noIC_deg8');
else
    load(filename);
    disp(['Processing ', filename]);
end

if nargin < 3
    xs = [0.35; 0; 0; 0.7];
end

if nargin < 4
    mycolor = [0 0 1];
end

PlotMySol_freeIC;

%% Trajecotry
axes(h_axis);
hold on;
axis equal;
xlim([-0.1,2.35]);
ylim([-0.1,0.6]);
plot(x_hist, y_hist,'b-','LineWidth', 2, 'color', mycolor);

%% Ground
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

%% Key frames
idx = [1; find( diff(mode_hist) == -1 )+1];
idx2 = [];
yR_lo = params.yR_lo;
for i = 1 : length(idx)-1
    tmp_y = y_hist( idx(i)+10: idx(i+1)-10 );
    [~, I] = min(abs(tmp_y - yR_lo));
    idx2 = [idx2; I+idx(i)+9];
end

if (current_mode == 0)
    idx = [idx; idx2];
else
    idx = [idx; idx2; length(t_hist)];
end

for i = 1 : length(idx)
    plotframe( gca, x_hist, y_hist, origin_hist, idx(i) );
end

xticks(unique(origin_hist));
