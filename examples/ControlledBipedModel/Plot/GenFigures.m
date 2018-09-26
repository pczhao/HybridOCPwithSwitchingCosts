% Generate figures
clear;
close all;

h_fig1 = figure(3);
subplot(4,1,1);
plotdata_freeIC(gca, 'Result_T3_freeIC_deg4', [ 0.35; 0; 0; 0.7 ], [.988,.553,.349]);

subplot(4,1,2);
plotdata_freeIC(gca, 'Result_T3_freeIC_deg6', [ 0.35; 0; 0; 0.73975 ], [.890,.290,.2]);

subplot(4,1,3);
plotdata_freeIC(gca, 'Result_T3_freeIC_deg8', [ 0.35; 0; 0; 0.72019 ], [.702,0,0]);

subplot(4,1,4);
plotdata_gpops( gca, 'Result_T3_freeIC_gpops' );

set(h_fig1,'PaperSize',[6, 7.4]);
h_fig1.PaperUnits = 'inches';
h_fig1.PaperPosition = [-0.4, -0.6, 6.8 8.4];
print(h_fig1,'Example1_freeIC','-dpdf');

% ============================================================================

h_fig2 = figure(4);
subplot(4,1,1);
plotdata_freeIC(gca, 'Result_T6_freeIC_deg4', [ 0.35; 0; 0; 0.7 ], [.988,.553,.349]);

subplot(4,1,2);
plotdata_freeIC(gca, 'Result_T6_freeIC_deg6', [ 0.35; 0; 0; 0.7 ], [.890,.290,.2]);

subplot(4,1,3);
plotdata_freeIC(gca, 'Result_T6_freeIC_deg8', [ 0.35; 0; 0; 0.7 ], [.702,0,0]);

subplot(4,1,4);
plotdata_gpops( gca, 'Result_T6_freeIC_gpops' );

set(h_fig2,'PaperSize',[6, 7.4]);
h_fig2.PaperUnits = 'inches';
h_fig2.PaperPosition = [-0.4, -0.6, 6.8 8.4];
print(h_fig2,'Example2_freeIC','-dpdf');
