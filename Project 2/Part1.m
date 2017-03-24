%Project #2: Part 1 - ME 303
%Sebastien Blanchet, Timothy Wulff

%Intialize script
close all
clear variables
clc

%Create symbolic variables for computation of integrals
syms n t x S L K T_L T_0 alpha ;

%Set integral functions, view equations in report
g(x)= ((S*((L^2)-(x^2))/(2*K)))+T_L;
p(x)= (T_0-g(x))*cos((n*pi*x)/(2*L));
q(x)= ((cos((n*pi*x)/(2*L)))^2);

%Perform symbolic inegration and display output to aid calculations
P(x)=int(p,0,L);
display(P);
Q(x)=int(q,0,L);
display(Q);

%Clear symbolic variables
clear variables

%Define constants
L = 0.15;               % length [m]
D_n = 0.01125;          % diameter [m]
alpha = 1.17e-4;        % thermal diffusivity [m^2/s]
P = 8;                  % input power [W]
row = 8933;             % density [kg/m^3]
C_p = 385;              % specific heat capacity [J/kg*deg.C]
K = 401;                % thermal conductivity [W/m]
T_0 = 18;               % outside temp. [deg.C]
T_L = 25;               % final temp. [deg.C]
t_end = 1000;           % end time [s]

%Constant Calculations
A = (pi/4)*(D_n^2);     % bar area [m^2]
S = P/(L*A);            % source strenght [W/m^3]
Q = S/(row*C_p);        % heat source term [deg.C/s]

%Calculate terms Dn terms for n=1,3,5 (first three non zero)
n=[1 3 5];

%Calculate D_n and lambda
lambda=(n*pi)/(2*L);
D_out= (((4*T_0)./(pi*n))-((4*T_L)./(pi*n))-((16*(L^2)*S)./(K*((pi.*n).^3)))).*(sin((n*pi)/(2)));


%Return Matrix of n, D_n , lambda and alpha*lambda^2 for Table 1
Output=[n;D_out;lambda;alpha*(lambda.^2)];

%Display result
format shortg
display(Output);

