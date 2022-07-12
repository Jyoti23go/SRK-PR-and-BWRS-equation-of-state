%         Benedict_Webb_Rubin_Starling equation    %
%%%To plot the variation in the value of compressibility factor with the Pressure(MPa) 
%......FOR METHANE.......%%%%%
z1=zeros(1,25);  %consider a zeros matrix of size 1*25 ,to store the value of z in it.
p1=zeros(1,25); %consider a zeros matrix of size 1*25 ,to store the value of P in it.
for i=1:25   %starting for loop
for p=1:25  %%input pressure in MPa
T=509.4;  % temperature IN Rankine %AT 283 KELVIN 
Tc=342.99;%% critical temperature in Rankine for methane
R=10.731573;   %%gas constant in ft^3.psia/R0.lb.mole
% all constant parameters for methane %%%
B0=0.723251;
A0=7520.29;
C0=2.71092*10^(8);
gama=1.48640;
b=0.925404;
a=2574.89;
alfa=0.468828;
c=4.37222*10^(8);
D0=1.07737*10^(10);
d=4.74891*10^(4);
E0=3.01122*10^(10);
P=145.03*p;         %%converting from MPa to psia to solve equation
%equation is:
% P=rho*R*T +(B0*R*T-A0-(C0/T^2)+(D0/T^3)-(E0/T^4))*(rho^2) +(b*R*T
% -a-(d/T))*(rho^3) +alfa*(a+(d/T))*(rho^6) +
% (c*rho^3/T^2)*(1+gama*(rho^2))*exp(-gama*(rho^2));
%%%%%to find the Taylor series of exponential term
%syms rho
%gama=1.48640;
%x= -(rho^2)*gama;
%T = taylor(exp(x))
%result for above is:
%(863041*rho^4)/781250 - (929*rho^2)/625 + 1
%finding coefficient for the rho variables:
a1=B0*R*T-A0-(C0/T^2)+(D0/T^3)-(E0/T^4);
a2=b*R*T-a-(d/T);
a3=alfa*(a+(d/T));
c1=(1.1*c*gama)/(T^2);
c2=0;
c3=(c/(T^2))*(1.1-1.48*gama);
c4=a3;
c5=(c/(T^2))*(gama-1.48);
c6=0;
c7=a2+(c/(T^2));
c8=a1;
c9=R*T;
c10=-P;
rho=roots([c1 c2 c3 c4 c5 c6 c7 c8 c9 c10]);
% Ignor all imaginary roots
rho=rho(imag(rho)==0);
%find z for the above rho
z=P/(rho*R*T);
p1(1,i)=p;   %storing value of 'P'
z1(1,i)=z;
i=i+1;
end
end
disp(z1);
plot(p1,z1,'-o');
axis([0 25 0.7 1]);
xlabel("Pressure in MPa");
ylabel("compressibility factor Z");
title('BWRS equation');