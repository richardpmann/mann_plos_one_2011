close all
marksize = 5;

load ./results/HUGEBENCHMARK_withBgen


%Make timestep figures


hold off
subplot(3,2, 1)
errorbar([1:10]-0.1, mean(meanR'), mean(stdR'), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(meanR2'), mean(stdR2'), 'ob', 'MarkerSize', marksize)
plot([0 11], 4*[1, 1], '--k')
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [3.6 4.9])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Number of Recorded Time-steps')
ylabel('Interaction Radius')
%legend('Initial 10 Timesteps', 'Steady State')
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'A')

hold off
subplot(3,2, 2)
errorbar(1:10, mean(meanB'), mean(stdB'), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(meanB2'), mean(stdB2'), 'ob', 'MarkerSize', marksize)
plot([0 11], 0.1*[1, 1], '--k')
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [0 0.22])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Number of Recorded Time-steps')
ylabel('Alignment')
legend('Initial', 'Steady State')
legend boxoff
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'B')


hold off
subplot(3,2, 3)
errorbar(1:10, mean(meanC'), mean(stdC'), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(meanC2'), mean(stdC2'), 'ob', 'MarkerSize', marksize)
plot([0 11], 1*[1, 1], '--k')
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [0.9 1.1])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Number of Recorded Time-steps')
ylabel('Attraction')
%legend('Initial 10 Timesteps', 'Steady State')
%text(0.1, 1.09, 'C')
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'C')

hold off
subplot(3,2, 4)
errorbar(1:10, mean(meanE'/pi), mean(stdE'/pi), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(meanE2'/pi), mean(stdE2'/pi), 'ob', 'MarkerSize', marksize)
plot([0 11], 0.0486*[1, 1], '--k')
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [0.04 0.07])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Number of Recorded Time-steps')
ylabel('Angular Noise / \pi')
%legend('Initial 10 Timesteps', 'Steady State')
%text(0.1, 0.078, 'D')
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'D')

hold off
subplot(3,2, 5)
errorbar(1:10, mean(meanBA'/pi), mean(stdBA'/pi), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(meanBA2'/pi), mean(stdBA2'/pi), 'ob', 'MarkerSize', marksize)
plot([0 11], (1/6)*[1, 1], '--k')
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [0.12 0.2])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Number of Recorded Time-steps')
ylabel('Blind Angle / \pi')
%legend('Initial 10 Timesteps', 'Steady State')
%text(0.1, 0.1945, 'E')
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'E')

hold off
subplot(3,2, 6)
errorbar(1:10, mean(Pentropy'/log(2)), std(Pentropy'/log(2)), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(Pentropy2'/log(2)), mean(Pentropy2'/log(2)), 'ob', 'MarkerSize', marksize)
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [0 36])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Number of Recorded Time-steps')
ylabel('Entropy / Bits')
%legend('Initial 10 Timesteps', 'Steady State')
%text(0.1, 38.3, 'F')
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(14/15)*diff(ylim), 'F')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.83 9*3/4]);

print -dpng -r300 time_multipanel.png
print -depsc2 time_multipanel.eps


%Make noise figures

close all

hold off
subplot(3, 2, 1)
errorbar(E_test/pi, mean(meanR_E'), mean(stdR_E'), '*k', 'MarkerSize', marksize)
hold on
plot([-.1,1.1], 4*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [3.5 4.8])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Simulation Angular Noise / \pi')
ylabel('Interaction Radius')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'A')


hold off
subplot(3, 2, 2)
errorbar(E_test/pi, mean(meanB_E'), mean(stdB_E'), '*k', 'MarkerSize', marksize)
hold on
plot([-.1 1.1], 0.1*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 0.5])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Simulation Angular Noise / \pi')
ylabel('Alignment')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'B')



hold off
subplot(3, 2, 3)
errorbar(E_test/pi, mean(meanC_E'), mean(stdC_E'), '*k', 'MarkerSize', marksize)
hold on
plot([-.1 1.1], 1*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0, 2])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%xlabel('Simulation Angular Noise / \pi')
ylabel('Attraction')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'C')


hold off
subplot(3, 2, 4)
errorbar(E_test/pi, mean(meanBA_E'/pi), mean(stdBA_E'/pi), '*k', 'MarkerSize', marksize)
hold on
plot([-.1 1.1], (1/6)*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0.1 0.2])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Simulation Angular Noise / \pi')
ylabel('Blind Angle / \pi')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'D')


% hold off
% errorbar(E_test/pi, mean(meanE_E'/pi), mean(stdE_E'/pi), '*k')
% hold on
% set(gca, 'XLim', [-.1, 1.1])
% set(gca, 'FontSize', 14)
% xlabel('Simulation Angular Noise / \pi')
% ylabel('Angular Noise / \pi')
% print -dpng -r300 postE_noise.png

hold off
subplot(3, 2, 5)
errorbar(E_test/pi, mean(Pentropy_E'/log(2)), std(Pentropy_E'/log(2)), '*k', 'MarkerSize', marksize)
hold on
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 25])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Simulation Angular Noise / \pi')
ylabel('Entropy / Bits')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'E')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.83 9*3/4]);

print -dpng -r300 noise_multipanel.png
print -depsc2 noise_multipanel.eps


%Make Q figures
close all

load ./results/HUGEBENCHMARK_Qfree *Q*
hold off
subplot(4, 2, 1)
errorbar(Q-0.01, mean(meanR_Q'), mean(stdR_Q'), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(meanR_Qfree'), mean(stdR_Qfree'), '^r', 'MarkerSize', marksize)
plot([-0.1 1.1], 4*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [3.5 4.9])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
legend('Fixed', 'Variable')
legend boxoff
%xlabel('Simulation Update Rate / Time-step')
ylabel('Interaction Radius')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'A')

hold off
subplot(4, 2, 2)
errorbar(Q-0.01, mean(meanB_Q'), mean(stdB_Q'), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(meanB_Qfree'), mean(stdB_Qfree'), '^r', 'MarkerSize', marksize)
plot([-.1 1.1], 0.1*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 0.2])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%legend('Fixed', 'Variable')
%xlabel('Simulation Update Rate / Time-step')
ylabel('Alignment')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'B')


hold off
subplot(4, 2, 3)
errorbar(Q-0.01, mean(meanC_Q'), mean(stdC_Q'), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(meanC_Qfree'), mean(stdC_Qfree'), '^r', 'MarkerSize', marksize)
plot([-.1 1.1], 1*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 1.2])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%legend('Update rate = 1', 'Variable update rate')
%xlabel('Simulation Update Rate / Time-step')
ylabel('Attraction')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'C')

hold off
subplot(4, 2, 4)
errorbar(Q-0.01, mean(meanE_Q'/pi), mean(stdE_Q'/pi), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(meanE_Qfree'/pi), mean(stdE_Qfree'/pi), '^r', 'MarkerSize', marksize)
plot([-.1 1.1], 0.05*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0.04 0.12])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%legend('Update rate = 1', 'Variable update rate')
%xlabel('Simulation Update Rate / Time-step')
ylabel('Angular Noise / \pi')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'D')

hold off
subplot(4, 2, 5)
errorbar(Q-0.01, mean(meanBA_Q'/pi), mean(stdBA_Q'/pi), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(meanBA_Qfree'/pi), mean(stdBA_Qfree'/pi), '^r', 'MarkerSize', marksize)
plot([-0.1 1.1], (1/6)*[1, 1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 0.3])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%legend('Update rate = 1', 'Variable update rate')
%xlabel('Simulation Update Rate / Time-step')
ylabel('Blind Angle / \pi')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'E')

hold off
subplot(4, 2, 6)
errorbar(Q-0.01, mean(Pentropy_Q'/log(2)), std(Pentropy_Q'/log(2)), '*k', 'MarkerSize', marksize)
hold on
errorbar(Q+0.01, mean(Pentropy_Qfree'/log(2)), std(Pentropy_Qfree'/log(2)), '^r', 'MarkerSize', marksize)
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 15])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
%legend('Update rate = 1', 'Variable update rate')
xlabel('Simulation Update Rate / Time-step')
ylabel('Entropy / Bits')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'F')


hold off
subplot(4, 2, 7)
errorbar(Q, mean(meanQ_Qfree, 2), mean(stdQ_Qfree, 2), '^r', 'MarkerSize', marksize)
hold on
plot([-0.1 1.1], [-0.1, 1.1], '--k')
set(gca, 'XLim', [-.1 1.1])
set(gca, 'YLim', [0 1.1])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Simulation Update Rate / Time-step')
ylabel('Update Rate')
ylim = get(gca, 'YLim');
text(-0.09, ylim(1)+(14/15)*diff(ylim), 'G')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.83 9]);

print -dpng -r300 Q_multipanel.png
print -depsc2 Q_multipanel.eps

%Make BF figures (requires loading other data files)

%Geo v Topo

close all


load ./results/HUGEBENCHMARK_withBgen BF BF2
hold off
subplot(1,2,1)
errorbar(1:10, mean(BF'/log(2)), std(BF'/log(2)), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(BF2'/log(2)), std(BF2'/log(2)), 'ob', 'MarkerSize', marksize)
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [-200 400])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Number of Recorded Time-steps')
ylabel('Log_2 Bayes Factor / Bits')

load ./results/HUGEBENCHMARK_withBtopogen BF BF2
errorbar(1:10, mean(BF'/log(2)), std(BF'/log(2)), '^r', 'MarkerSize', marksize)
errorbar([1:10]+0.1, mean(BF2'/log(2)), std(BF2'/log(2)), '.g', 'MarkerSize', marksize)
legend('Initial, G', 'Steady State, G', 'Initial, T', 'Steady State, T', 'Location', 'NorthWest')
legend boxoff
axis square
ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(0.5/15)*diff(ylim), 'A')
%lh1 = legend('Initial 10 timesteps, geometric model', 'Steady State, geometric model', 'Initial 10 timesteps, topological model', 'Steady State, topological model', 'Location', 'NorthWest')
%set(lh1, 'box', 'off')

load ./results/HUGEBENCHMARK_withBgen BFB BFB2

hold off
subplot(1,2,2)

errorbar(1:10, mean(BFB'/log(2)), std(BFB'/log(2)), '*k', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(BFB2'/log(2)), std(BFB2'/log(2)), 'ob', 'MarkerSize', marksize)
set(gca, 'XLim', [0 11])
set(gca, 'YLim', [-20 100])
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
xlabel('Number of Recorded Time-steps')
ylabel('Log_2 Bayes Factor / Bits')



load ./results/HUGEBENCHMARK_withoutBgen BFB BFB2
errorbar(1:10, mean(BFB'/log(2)), std(BFB'/log(2)), '^r', 'MarkerSize', marksize)
hold on
errorbar([1:10]+0.1, mean(BFB2'/log(2)), std(BFB2'/log(2)), '.g', 'MarkerSize', marksize)

legend('Initial, +A', 'Steady State, +A', 'Initial, -A', 'Steady State, -A', 'Location', 'NorthWest')
legend boxoff
axis square
%pos2 = get(lh2, 'Position');
%set(lh2, 'box', 'off', 'Position', pos2 - [0.05, 0, 0, 0])

ylim = get(gca, 'YLim');
text(0.1, ylim(1)+(0.5/15)*diff(ylim), 'B')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.83 9/4]);



print -dpng -r300 bayes_factor_multipanel.png
print -depsc2 bayes_factor_multipanel.eps

close all

figure 
subplot(1,2,1)
sppABC(10, 25, 4, 1, 1, 1, 0.1, 1, 0.1*pi, pi/6, 1, 1);
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
text(0.1, 0.3, 'A')
axis square
subplot(1,2,2)
sppABC(10, 25, 4, 1, 50, 1, 0.1, 1, 0.1*pi, pi/6, 1, 1);
set(gca, 'FontSize', 10, 'FontName', 'Ariel')
text(0.1, 0.3, 'B')
axis square
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.83 3.7]);

print -dpng -r300 torus_multipanel.png
print -depsc2 torus_multipanel.eps




