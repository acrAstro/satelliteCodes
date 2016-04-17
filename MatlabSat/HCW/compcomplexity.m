clear all; close all; clc;

t = [1 10 100 1000 10000];

S1 = [0.0161 0.0433 0.3871 3.8009 37.8485];
S2 = [0.2124 0.4651 3.3159 32.1149 319.3837];

figure(1)
hold on
grid on
semilogx(t,S1,'k','LineWidth',2)
semilogx(t,S2,'r','LineWidth',2)
title1 = title('Run Time vs Number of Orbits');
xl = xlabel('Number of Orbits, $N$');
yl = ylabel('CPU Time, $s$');
leg1 = legend('STM','DNS','Location','Best');
set([leg1 title1 xl yl],'interpreter','latex','fontsize',12)
axis tight