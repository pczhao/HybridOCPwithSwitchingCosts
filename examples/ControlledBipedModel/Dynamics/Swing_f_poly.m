function [out] = Swing_f_poly( in1, params )

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
  [mydq1;
  0.E-323+0.1E1.*mydq2.^2.*myq1+(1/2).*g.*((-2)+myq2.^2);
  [mydq2];
  l0.^(-3).*(0.1E1.*g.*myq1.^2.*myq2+l0.*(0.2E1.*mydq1.* ...
  mydq2.*myq1+(-0.3E1).*g.*myq1.*myq2)+l0.^2.*((-0.4E1).*mydq1.* ...
  mydq2+g.*myq2.*(0.3E1+(-0.166667E0).*myq2.^2)))];