%%%%%%%%%%%........Soave-Redlich-Kwong...(SRK).....Equation......%%%%%%%%
%%%To plot the variation in the value of compressibility factor with the Pressure(MPa) 
%......FOR HYDROGEN.......%%%%%
Tc=33.2;    %critical temperature in kelvin for hydrogen
Pc=1.3;     %critical pressure in 'MPa' for hydrogen
T=283;      % considered temperature T=283K
omega=-0.218; %accentric factor for hydrogen
z=zeros(1,50);%consider a zeros matrix of size 1*50,to store the value of z in it.
p=zeros(1,50);%consider a zeros matrix of size 1*50,to store the value of pressure in it.
for j=1:50%started a for loop to store all values
    
for P=1:1:50   %started another for loop to find values of 'Z' at different Pressure
Tr=T/Tc;        %%Reduce temperature
Pr=P/Pc;        %%Reduce Pressure
m=0.48508+1.55171*omega-0.15613*omega^2;  %%%fixed for SRK equation
alfa=(1+m*(1-sqrt(Tr)))^2;        %%value of alfa
A=0.427*alfa*Pr/Tr;              %%value of parameter 'A'
B=0.08664*Pr/Tr;                 %%value of parameter 'B'
%%%   z^3-z^2+z(A-B-B^2)-A*B==0       %%cubic equation for SRK
a=1;    %coefficient for z^3
b=-1;   %coefficient for z^2
c=A-B-B^2;  %coefficient for z
d=-A*B;     %constant term
%%%%P=R*T/(V-b)+a*alfa/(V(V+b));%%%%FOR SRK
z1=roots([a,b,c,d]);    %%to find the roots of cubic equation
z1=z1(imag(z1)==0);     %%Neglecting the imaginary roots,consider only real root(out of 3 roots of cubic equation)
p(1,j)=P;   %storing value of 'P'
z(1,j)=z1;  %storing value of 'z'
j=j+1;      %increasing value of j by 1;
end
end
disp(z); %to display the compressibility factor values
n=size(z);  %%to find the size of matrix of z 
disp('size of z');  %%to display the size 
disp(n);                %%to display the size 
disp(p);                %%to display value of pressure
plot(p,z,'-o'); %%to plot the graph for p v/s z
axis([0 50 0.7 1]);       %%to give axis of graph
grid on 
title("For SRK ")%% grid on
xlabel('Pressure in MPa') %%labelling x axis
ylabel('value of Z')     %%labelling y axis
