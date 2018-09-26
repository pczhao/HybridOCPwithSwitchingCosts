% Generate figures
clear;
close all;

h_fig1 = figure(3);
subplot(4,1,1);
plotdata_freeIC(gca, 'Data/Result_T3_freeIC_deg4', [ 0.35; 0; 0; 0.7 ], [.988,.553,.349]);

subplot(4,1,2);
plotdata_freeIC(gca, 'Data/Result_T3_freeIC_deg6', [ 0.35; 0; 0; 0.7432 ], [.890,.290,.2]);

subplot(4,1,3);
plotdata_freeIC(gca, 'Data/Result_T3_freeIC_deg8', [ 0.35; 0; 0; 0.72792 ], [.702,0,0]);

subplot(4,1,4);
plotdata_gpops( gca, 'Data/Result_T3_freeIC_gpops' );

suptitle('T=3');

% ============================================================================

h_fig2 = figure(4);
subplot(4,1,1);
plotdata_freeIC(gca, 'Data/Result_T6_freeIC_deg4', [ 0.35; 0; 0; 0.7 ], [.988,.553,.349]);

subplot(4,1,2);
plotdata_freeIC(gca, 'Data/Result_T6_freeIC_deg6', [ 0.35; 0; 0; 0.7 ], [.890,.290,.2]);

subplot(4,1,3);
plotdata_freeIC(gca, 'Data/Result_T6_freeIC_deg8', [ 0.35; 0; 0; 0.7 ], [.702,0,0]);

subplot(4,1,4);
plotdata_gpops( gca, 'Data/Result_T6_freeIC_gpops' );

suptitle('T=6');