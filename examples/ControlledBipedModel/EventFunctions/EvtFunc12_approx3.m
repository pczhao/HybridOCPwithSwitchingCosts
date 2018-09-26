function [value, isterminal, direction] = EvtFunc12_approx3( ~, x, params )
% {y <= yR} -> {y >= yR}

polysin = @(ang) ang - ang^3/6;
polycos = @(ang) 1 - ang^2/2;
% polysin = @(ang) ang - ang.^3/6 + ang.^5/120;
% polycos = @(ang) 1 - ang.^2/2 + ang.^4/24;

y = x(1) * polycos(x(3));
yR_hi = params.yR_hi;
yR_lo = params.yR_lo;

isterminal = [1;1];
direction = [1;0];
% value = y - yR_hi;
value = [ y - yR_lo; y];