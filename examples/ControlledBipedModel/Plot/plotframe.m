function plotframe(h_axis, x_hist, y_hist, origin_hist, idx)

x = x_hist( idx );
y = y_hist( idx );
O = origin_hist( idx );

% Draw leg
plot(h_axis, [x,O], [y,0], 'k-', 'LineWidth', 2);

% Draw mass
Mpos = [x,y];
th = 0:pi/50:2*pi;
xvec = Mpos(1) + 0.03*cos(th);
yvec = Mpos(2) + 0.03*sin(th);
patch(xvec, yvec, 1, 'FaceColor', [0.8,0.8,0.8], 'LineWidth',2);

