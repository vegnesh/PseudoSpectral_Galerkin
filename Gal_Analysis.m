clc
clear all 
close all
%% Galerikin Results Data Analysis

load('Galerkin_data.mat');
Gam = 1.4;
R = 1/Gam;
xc = linspace(0,2*pi,1000);
colorv = ['b','k','r','y','g','c','m'];
index5e3 = [1,21,41,61,81,101,106];
index1e3 = [1,101,201,301,401,501,541];
index5e4 = [1,201,401,601,801,1001,1010];
for ival = length(index5e3)
    Nf = 90;
    tempN = 2*Nf+1;
    rk = Qfsoln90_5e3(1:tempN,index5e3(ival));
    uk = Qfsoln90_5e3(tempN+1:2*tempN,index5e3(ival));
    pk = Qfsoln90_5e3(2*tempN+1:end,index5e3(ival));
    lw2 =2.5;
    figure(1);
    subplot(2,2,1)
    plot(xc/pi,1./(getfunval(rk,Nf,xc)),colorv(ival),'Linewidth',lw2)
    hold on 
    title('Density')
    xlabel('x')
    subplot(2,2,2)
    plot(xc/pi,getfunval(pk,Nf,xc),colorv(ival),'Linewidth',lw2)
    title('Pressure')
    xlabel('x')
    hold on 
    subplot(2,2,3)
    plot(xc/pi,getfunval(pk,Nf,xc).*getfunval(rk,Nf,xc)/R,colorv(ival),'Linewidth',lw2)
    title('Temperature')
    xlabel('x')
    hold on 
    subplot(2,2,4)
    plot(xc/pi,getfunval(uk,Nf,xc),colorv(ival),'Linewidth',lw2)
    hold on 
    title('Velocity')
    xlabel('x')
    
    Nf = 60;
    tempN = 2*Nf+1;
    rk = Qfsoln60_5e3(1:tempN,index5e3(ival));
    uk = Qfsoln60_5e3(tempN+1:2*tempN,index5e3(ival));
    pk = Qfsoln60_5e3(2*tempN+1:end,index5e3(ival));
    lw2 =2;
    figure(1);
    subplot(2,2,1)
    plot(xc/pi,1./(getfunval(rk,Nf,xc)),colorv(ival-5),'Linewidth',lw2)
    hold on 
    title('Density')
    xlabel('x')
    subplot(2,2,2)
    plot(xc/pi,getfunval(pk,Nf,xc),colorv(ival-5),'Linewidth',lw2)
    title('Pressure')
    xlabel('x')
    hold on 
    subplot(2,2,3)
    plot(xc/pi,getfunval(pk,Nf,xc).*getfunval(rk,Nf,xc)/R,colorv(ival-5),'Linewidth',lw2)
    title('Temperature')
    xlabel('x')
    hold on 
    subplot(2,2,4)
    plot(xc/pi,getfunval(uk,Nf,xc),colorv(ival-5),'Linewidth',lw2)
    hold on 
    title('Velocity')
    xlabel('x')
    legend('90 modes at t=0.5','60 modes at t=0.5','Location','South')
end    
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [11 8]);
set(gcf,'PaperPosition',[0.5,0.5,10,7]);
set(gcf,'PaperPositionMode','Manual');
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
print(gcf, '-dpdf', '-r150', 'GS_Temp.pdf');

for ival = 1:length(index5e3)
    Nf = 90;
    tempN = 2*Nf+1;
    rk = Qfsoln90_5e3(1:tempN,index5e3(ival));
    uk = Qfsoln90_5e3(tempN+1:2*tempN,index5e3(ival));
    pk = Qfsoln90_5e3(2*tempN+1:end,index5e3(ival));
    lw2 =2.5;
    figure(2);
    subplot(2,2,1)
    plot(xc/pi,1./(getfunval(rk,Nf,xc)),colorv(ival),'Linewidth',lw2)
    hold on 
    title('Density')
    xlabel('x')
    subplot(2,2,2)
    plot(xc/pi,getfunval(pk,Nf,xc),colorv(ival),'Linewidth',lw2)
    title('Pressure')
    xlabel('x')
    hold on 
    subplot(2,2,3)
    plot(xc/pi,getfunval(pk,Nf,xc).*getfunval(rk,Nf,xc)/R,colorv(ival),'Linewidth',lw2)
    title('Temperature')
    xlabel('x')
    hold on 
    subplot(2,2,4)
    plot(xc/pi,getfunval(uk,Nf,xc),colorv(ival),'Linewidth',lw2)
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
print(gcf, '-dpdf', '-r150', 'GS_Cumulative.pdf');
