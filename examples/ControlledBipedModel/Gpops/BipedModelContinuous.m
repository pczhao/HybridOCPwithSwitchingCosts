%-----------------------------------------%
% BEGIN: BipedModelContinuous.m %
%-----------------------------------------%
function phaseout = BipedModelContinuous(input)

nphases = input.auxdata.nphases;
params = input.auxdata.params;
T = input.auxdata.T;
phaseout = struct('dynamics',cell(1,nphases),'integrand',cell(1,nphases));

polysin = @(ang) ang - ang.^3/6;
polycos = @(ang) 1 - ang.^2/2;

% keyboard;
for iphase = 1 : nphases
    idx = mod( iphase+1, 2 ) + 1;
    t1 = input.phase(iphase).time;
    x1 = input.phase(iphase).state;
    u1 = input.phase(iphase).control;
    
    y = x1(:,1) .* polycos(x1(:,3));
    switch idx
        case 1          % Stance phase, y<=yR
            phaseout(iphase).dynamics = ( Swing_f_poly(x1', params) + Swing_g_poly(x1', params) * u1' )';
            phaseout(iphase).integrand = (u1) .^ 2;
            phaseout(iphase).path = y;
        case 2          % Stance phase, y>=yR
            phaseout(iphase).dynamics = ( Swing_f_poly(x1', params) + Swing_g_poly(x1', params) * u1' )';
            phaseout(iphase).integrand = (u1) .^ 2;
            phaseout(iphase).path = y;
        otherwise
            error('Wrong phase.');
    end
    
end

%---------------------------------%
% END: BipedModelContinuous.m %
%---------------------------------%