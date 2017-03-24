%Project #2: Part 5 - ME 303

%Sebastien Blanchet, Timothy Wulff
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
dx = 0.01;              % delta x [m]
dt=0.125;                   % delta t [s]

%Calculations
A = (pi/4)*(D_n^2);     % bar area [m^2]
S = P/(L*A);            % source strenght [W/m^3]
Q = S/(row*C_p);        % heat source term [deg.C/s]

%Analytical solution for T(x,t) for x at t=40s
t_n = 40;
x = 0:dx:L;

%Create for loop to evalute infinite sum numerically
for i=1:L/dx+1
    for n=1:2:1001
        D_n(n,i)=(((4*T_0)/(pi*n))-((4*T_L)/(pi*n))-((16*(L^2)*S)/(K*((pi*n)^3))))*(sin((n*pi)/(2)));
        phi_x(n,i)= cos((n*pi*x(i))/(2*L));
        tau_t(n,i)= exp((-alpha*(((n*pi)/(2*L))^2)*(t_n)));
        f_xt(n,i)=D_n(n,i)*phi_x(n,i)*tau_t(n,i);
    end
end

%Add steady state solution to summation of all f(x,t) values
T_ss1=((-S/(2*K))*(x.^2))+(T_L+((S*(L^2))/(2*K)));
T_an1=sum(f_xt)+(T_ss1);
Tn=GetTnum(dt);

%Analytical solution for T(x,t) for 0<t<1000 at x=0
t = 0:dt:t_end;
x_n = 0;

%Create for loop to evalute infinite sum numerically
for i=1:t_end/dt+1
    for n=1:2:1001
        D_n(n,i)=(((4*T_0)/(pi*n))-((4*T_L)/(pi*n))-((16*(L^2)*S)/(K*((pi*n)^3))))*(sin((n*pi)/(2)));
        phi_x(n,i)= cos((n*pi*x_n)/(2*L));
        tau_t(n,i)= exp((-alpha*(((n*pi)/(2*L))^2)*(t(i))));
        f_xt(n,i)=D_n(n,i)*phi_x(n,i)*tau_t(n,i);
    end
end

%Add steady state solution to summation of all f(x,t) values
T_ss2=((-S/(2*K))*(x_n^2))+(T_L+((S*(L^2))/(2*K)));
T_an2=sum(f_xt)+(T_ss2);


%Analytical vs numerical plot of T(x,t) at t=40s 
figure3=figure;
hold on
plot(x,T_an1);
plot(x,Tn(:,t_n/dt+1));
xlabel('L [m]');
ylabel('T [deg.C]');
title('Plot of T vs L at 40s');
legend('Analytical','Numerical');

%Analytical vs numerical plot of T(x,t) at steady state
figure4=figure;
hold on
plot(x,T_ss1);
plot(x,Tn(:,t_end/dt+1));
xlabel('L [m]');
ylabel('T [deg.C]');
title('Plot of T vs L at Steady State');
legend('Analytical','Numerical');

%Analytical vs numerical plot of T(x,t) at x=0
figure5=figure;
hold on
plot(t,T_an2);
plot(t,Tn(1,:));
xlabel('t [s]');
ylabel('T [deg.C]');
title('Plot of T vs t at x=0');
legend('Analytical','Numerical');
