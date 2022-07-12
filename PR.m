%%%%%%%......PengRobinson...(PR)...equation.......%%%%%%%
 
z2=zeros(1,50); %consider a zeros matrix of size 1*50,to store the value of z2 in it.
p1=zeros(1,50); %consider a zeros matrix of size 1*50 ,to store the value of pressure in it
for i=1:50       %started a for loop to store all values
for P1=1:1:50       %started another for loop to find values of 'Z' at different Pressure
Tr1=T/Tc;       %Reduce Temperature
Pr1=P1/Pc;          %Reduce Pressure
m1=0.37464+1.54226*omega-0.26992*omega^2;   %m1 for pr equation
alfa1=(1+m1*(1-sqrt(Tr1)))^2;               %alfa1 for pr equation
A1=0.45724*alfa*Pr1/Tr1;                    % 'A1' for pr equation
B1=0.07780*Pr1/Tr1;                         % 'B1' for pr equation
%%%   z^3-(1-B1)*z^2+(A1-3*B1^2-2*B1)*z-(A1*B1-B1^2-B1^3)==0;       %%cubic
%%%   equation for peng robinson equation
a1=1;       %coefficient for z^3
b1=B1-1;        %coefficient for z^2
c1=A1-3*(B1^2)-2*B1;        %coefficient for z
d1=-A1*B1+B1^2+B1^3;        %constant term
%%%%% P=R*T/(V-b)-a*alfa/(V*(V+b)+b*(V-b));
z3=roots([a1,b1,c1,d1]);        %%to find the roots of cubice equation
z3=z3(imag(z3)==0);             %%Neglecting the imaginary roots,consider only real root(out of 3 roots of cubic equation)
z2(1,i)=z3;          %storing value of 'z'
p1(1,i)=P1;              %storing value of 'Pressure'
i=i+1;                  %increasing value of i by 1
end
end
disp(z3);           %display the value for z3
m=size(z3);             %find the size of z3 matrix
disp('size of z3');         %display the size of z3 matrix
disp(m);                        %display the size of z3 matrix
disp(p1);                   %display the value for p1
plot(p1,z2,'-s');           %plot for p1 v/s z2
axis([0 50 0.7 1]);         %give the axis limit
xlabel('Pressure');
ylabel('Z-factor');
title('Peng Robinson Equation');