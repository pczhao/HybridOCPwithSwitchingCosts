%---------------------------------------%
% BEGIN: BipedModelEndpoint.m %
%---------------------------------------%
function output = BipedModelEndpoint(input)

nphases = input.auxdata.nphases;
params = input.auxdata.params;
l0 = input.auxdata.l0;
yR_lo = input.auxdata.yR_lo;
yR_hi = input.auxdata.yR_hi;
al = input.auxdata.alpha;
T = input.auxdata.T;

polysin = @(ang) ang - ang.^3/6;
polycos = @(ang) 1 - ang.^2/2;

objective = 0;

% Events
for iphase = 1 : nphases-1
    idx = mod( iphase+1, 2 ) + 1;
    
    switch idx
        case 1      % Stance phase, y<=yR_hi
            tf1 = input.phase(iphase).finaltime;
            xf1 = input.phase(iphase).finalstate;
            t02 = input.phase(iphase+1).initialtime;
            x02 = input.phase(iphase+1).initialstate;
            
            y = xf1(1) * polycos(xf1(3));
            G = y - yR_hi;                      % guard
            R = xf1;                            % reset map
            output.eventgroup(iphase).event = [ t02 - tf1, x02 - R, G ];
            
        case 2      % Stance phase, y>=yR_lo
            tf1 = input.phase(iphase).finaltime;
            xf1 = input.phase(iphase).finalstate;
            t02 = input.phase(iphase+1).initialtime;
            x02 = input.phase(iphase+1).initialstate;
            
            y = xf1(1) * polycos(xf1(3));
            G = y - yR_lo;                      % guard
            R = Reset_poly( xf1', params )';    % reset map
            output.eventgroup(iphase).event = [ t02 - tf1, x02 - R, G ];
            
            SwitchingCost = 10 * (xf1(1) * polysin(xf1(3)) + l0 * polysin(-al) - input.auxdata.d_des)^2;
            objective = objective + SwitchingCost;

    end
end

% Terminal condition
iphase = nphases;
tf1 = input.phase(iphase).finaltime;
output.eventgroup(iphase).event = [ tf1 - T ];

% Objective function
for i = 1 : nphases
    objective = objective + input.phase(i).integral;
end

output.objective = objective;

%---------------------------------%
% END: BipedModelEndpoint.m %
%---------------------------------%