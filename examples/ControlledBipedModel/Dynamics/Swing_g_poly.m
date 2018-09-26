function [out] = Swing_g_poly( in1, params )

if nargin < 2
    m = 1;
    g = 9.8;
    l0 = 1;
    error('No params!');
else
    m = params.m;
    g = params.g;
    l0 = params.l0;
end

myq1    = in1(1,:);
mydq1   = in1(2,:);
myq2    = in1(3,:);
mydq2   = in1(4,:);

out = ...
  [ 0
    1/m
    0
    0 ];
% out = repmat( out, 1, size( myq1, 2 ) );