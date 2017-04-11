clc
%clear all
close all

%% Galerkin

Nf = 50;

funct = @(x) exp(-(x-pi).^2);
dfunct= @(x) exp(-(x-pi).^2).*(-2*(x-pi));
A = getcoeff(funct,Nf,1000);

xc = linspace(0,2*pi,1000);
ya = funct(xc);
%plot(xc,dfunct(xc),'Linewidth',2);

Anew = [-Nf:Nf]';
Anew = 1i*Anew;
Anew = Anew.*A;
%hold on 
yf = getfunval(Anew,Nf,xc);
%plot(xc,yf,'r--','Linewidth',2);

%% U*DU

Ai = getcoeff(funct,Nf,1000);
kvec = [-Nf:Nf]';
Asol = ode45(@(t,A) -conv(A,A.*1i.*kvec,'same'),[0:0.1:1],Ai);

Acoeff = Asol.y;

for ival = length(Acoeff(1,:))
    plot(xc,getfunval(Acoeff(:,ival),Nf,xc),'Linewidth',2);
    hold on 
     plot(xc2,soln(:,ival),'r','Linewidth',2)
    pause(0.1);
    hold off
    

end

