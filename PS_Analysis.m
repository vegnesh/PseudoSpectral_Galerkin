clc
clear all 
close all

%% Data Analysis

load('PS_data.mat');
Gam = 1.4;
R = 1/Gam;
index1e3 = [1,101,201,301,401,501,530];
index1e4 = [1,1001,2001,3001,4001,5001,5301];

for ival = length(index1e3)
TwoN = length(Soln500_1e3(:,1))/3;
LRho = Soln500_1e3(1:TwoN,index1e3(ival));
LRU = Soln500_1e3(TwoN+1:2*TwoN,index1e3(ival));
LRE = Soln500_1e3(2*TwoN+1:end,index1e3(ival));

LU = LRU./LRho;
LP = (Gam-1)*(LRE - 0.5*LRU.*LU);
LT = LP./LRho/R;
xc = xc500_1e3;
figure(1)
lw1 =1;
subplot(2,2,1)
plot(xc,LRho,'r','Linewidth',lw1)
hold on 
title('Density')
subplot(2,2,2)
plot(xc,LP,'r','Linewidth',lw1)
title('Pressure')
hold on 
subplot(2,2,3)
plot(xc,LT,'r','Linewidth',lw1)
title('Temp')
hold on 
subplot(2,2,4)
plot(xc,LU,'r','Linewidth',lw1)
hold on 
title('Velocity')


TwoN = length(Soln1000_1e3(:,1))/3;
LRho = Soln1000_1e3(1:TwoN,index1e3(ival));
LRU = Soln1000_1e3(TwoN+1:2*TwoN,index1e3(ival));
LRE = Soln1000_1e3(2*TwoN+1:end,index1e3(ival));

LU = LRU./LRho;
LP = (Gam-1)*(LRE - 0.5*LRU.*LU);
LT = LP./LRho/R;
lw2 =1;
xc = xc1000_1e3;
figure(1)
subplot(2,2,1)
plot(xc,LRho,'b','Linewidth',lw2)
hold on 
title('Density')
xlabel('x')
subplot(2,2,2)
plot(xc,LP,'b','Linewidth',lw2)
title('Pressure')
xlabel('x')
hold on 
subplot(2,2,3)
plot(xc,LT,'b','Linewidth',lw2)
title('Temperature')
xlabel('x')
hold on 
subplot(2,2,4)
plot(xc,LU,'b','Linewidth',lw2)
hold on 
title('Velocity')
xlabel('x')
legend('500 points , dt = 1e-3','1000 points , dt = 1e-3','Location','South')
end
%legend('t=0','t=0.1','t=0.2','t=0.3','t=0.4','t=0.5','t=0.53','Location','BestOutside')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [11 8]);
set(gcf,'PaperPosition',[0.5,0.5,10,7]);
set(gcf,'PaperPositionMode','Manual');
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
print(gcf, '-dpdf', '-r150', 'ConvergencePS2.pdf');

colorv = ['b','k','r','y','g','c','m'];
for ival =1: length(index1e3)
TwoN = length(Soln500_1e3(:,1))/3;
LRho = Soln500_1e3(1:TwoN,index1e3(ival));
LRU = Soln500_1e3(TwoN+1:2*TwoN,index1e3(ival));
LRE = Soln500_1e3(2*TwoN+1:end,index1e3(ival));

LU = LRU./LRho;
LP = (Gam-1)*(LRE - 0.5*LRU.*LU);
LT = LP./LRho/R;
xc = xc500_1e3;

figure(2)
subplot(2,2,1)
plot(xc,LRho,colorv(ival),'Linewidth',lw2)
hold on 
title('Density')
xlabel('x')
subplot(2,2,2)
plot(xc,LP,colorv(ival),'Linewidth',lw2)
title('Pressure')
xlabel('x')
hold on 
subplot(2,2,3)
plot(xc,LT,colorv(ival),'Linewidth',lw2)
title('Temperature')
xlabel('x')
hold on 
subplot(2,2,4)
plot(xc,LU,colorv(ival),'Linewidth',lw2)
hold on 
title('Velocity')
xlabel('x')





end
legend('t=0','t=0.1','t=0.2','t=0.3','t=0.4','t=0.5','t=0.53','Location','BestOutside')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [11 8]);
set(gcf,'PaperPosition',[0.5,0.5,10,7]);
set(gcf,'PaperPositionMode','Manual');
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
print(gcf, '-dpdf', '-r150', 'PS_Cumulative.pdf');