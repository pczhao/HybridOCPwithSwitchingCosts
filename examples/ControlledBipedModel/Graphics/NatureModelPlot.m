classdef NatureModelPlot < handle
    
    properties (SetAccess = public)
        mass_hist;
        time_hist;
        O_hist;
        current_q;      % [ x, y, l, theta ], where (x,y) is the position of mass
        origin;       % x-position of the origin O
        phase;
        Mpos;
        Opos;
        params;
        mycolor;
        mythickness;
        
        flight_flag;
    end
    
    properties (SetAccess = private)
        h_2d;
        h_M;
        h_O;
        h_link;
        h_link2;
        h_gnd;
        h_traj;
    end
    
    methods
        % Constructor
        function P = NatureModelPlot(phase, x0, params, mycolor, mythickness)
            P.phase = phase;
            P.params = params;
            P.flight_flag = 0;
            P.Opos = [ 0, 0 ];
            if nargin > 3
                P.mycolor = mycolor;
                P.mythickness = mythickness;
            else
                P.mycolor = [];
            end
            switch P.phase
                case {1,2}  % Stance, y<= yR
                    ll = x0(1);
                    tt = x0(3);
                    P.current_q(1) = P.Opos(1) + ll * sin(tt);
                    P.current_q(2) = P.Opos(2) + ll * cos(tt);
                    P.current_q(3) = ll;
                    P.current_q(4) = tt;
                    P.O_hist = 0;
                case 3      % Flight
                    error('Should not be here.');
%                     P.current_q(1) = x0(1);
%                     P.current_q(2) = x0(3);
%                     P.current_q(3) = params.l0;
%                     P.current_q(4) = params.alpha;
                otherwise
                    error('Phase index error.');
            end
            
            figure(1);
%             axis equal;
            xlim([-0.5,5]);
            ylim([-0.1,1]);
            box on;
            P.h_2d = gca;
            hold( P.h_2d, 'on' );
            P.UpdatePos;
            P.h_M = plot(P.h_2d, ...
                         P.Mpos(1), P.Mpos(2), 'ro', 'MarkerSize', 10);
            P.h_O = plot(P.h_2d, ...
                         P.O_hist, P.O_hist * 0, 'ko', 'MarkerSize', 10);
            if isempty(P.mycolor)
                P.h_traj = plot(P.h_2d, [P.Mpos(1)], [P.Mpos(2)], '-r' );
            else
                P.h_traj = plot(P.h_2d, [P.Mpos(1)], [P.Mpos(2)], 'color', P.mycolor, 'LineWidth', P.mythickness);
            end
            P.h_link2 = plot( P.h_2d, ...
                              [P.Mpos(1), P.Mpos(1)-params.l0 * sin(params.alpha)], ...
                              [P.Mpos(2), P.Mpos(2)-params.l0 * cos(params.alpha)], ...
                              'Color', [0.4,0.4,0.4], 'LineWidth', 4 );
            P.h_link = plot(P.h_2d, ...
                            [P.Mpos(1), P.Opos(1)], ...
                            [P.Mpos(2), P.Opos(2)], ...
                            '-b', 'LineWidth', 4 );
            P.h_gnd = plot(P.h_2d, ...
                           [-2, 13], [0,0], '-k', 'LineWidth', 2);
            
            hold( P.h_2d, 'off' );
            title( P.h_2d, 0 );
        end
        
        function UpdatePlot( P, tval, xval )
            switch P.phase
                case {1,2}  % Stance
                    ll = xval(1);
                    tt = xval(3);
                    P.current_q(1) = P.Opos(1) + ll * sin(tt);
                    P.current_q(2) = P.Opos(2) + ll * cos(tt);
                    P.current_q(3) = ll;
                    P.current_q(4) = tt;
                case 3      % Flight
                    error('Should not be here.');
                    P.current_q(1) = xval(1);
                    P.current_q(2) = xval(3);
                    P.current_q(3) = P.params.l0;
                    P.current_q(4) = P.params.alpha;
                otherwise
                    error('Phase index error.');
            end
            P.UpdatePos;
            set( P.h_M, ...
                 'XData', P.Mpos(1), ...
                 'YData', P.Mpos(2) );
            set( P.h_O, ...
                 'XData', P.O_hist, ...
                 'YData', P.O_hist * 0 );
            if (P.flight_flag == 1) && (xval(1) > P.params.l0)
                set( P.h_link, ...
                 'XData', [P.Mpos(1), P.Mpos(1)-P.params.l0 * sin(P.params.alpha)], ...
                 'YData', [P.Mpos(2), P.Mpos(2)-P.params.l0 * cos(P.params.alpha)] );
            else
                set( P.h_link, ...
                 'XData', [P.Mpos(1), P.Opos(1)], ...
                 'YData', [P.Mpos(2), P.Opos(2)] );
            end
            set( P.h_link2, ...
                 'XData', [P.Mpos(1), P.Mpos(1)-P.params.l0 * sin(P.params.alpha)], ...
                 'YData', [P.Mpos(2), P.Mpos(2)-P.params.l0 * cos(P.params.alpha)] );
            set( P.h_traj, ...
                 'XData', P.mass_hist(:,1), ...
                 'YData', P.mass_hist(:,2) );
            title( P.h_2d, tval );
            
%             set( P.h_M, 'Marker', 'none' );
%             set( P.h_O, 'Marker', 'none' );
%             set( P.h_link, 'LineStyle', 'none' );
%             set( P.h_link2, 'LineStyle', 'none' );
            
            drawnow;
        end
            
        function UpdatePos( P )
%             ll = P.current_q(3);
%             tt = P.current_q(4);
            P.Mpos = P.current_q(1:2);
%             P.Opos = P.Mpos + [ll*sin(tt), -ll*cos(tt)];
            P.mass_hist = [P.mass_hist; P.Mpos(1:2)];
        end
        
        function Visualize( P, t_hist, x_hist, phase )
            P.phase = phase;
            for i = 2 : length(t_hist)
                tic;
                tval = t_hist(i);
                xval = x_hist(i,:);
                P.UpdatePlot( tval, xval );
                tstep = 0.5 * (t_hist(i) - t_hist(i-1));
                t = toc;
                if (t < tstep)
                    pause(tstep - t);
                end
            end
            ll = xval(1);
            tt = xval(3);
            if phase == 2   % 2 -> 1
                P.Opos = P.Opos + [ ll * sin(tt) + P.params.l0 * sin(-P.params.alpha), 0 ];
                P.O_hist = [ P.O_hist; P.Opos(1) ];
            end
        end
    end
end
            
            