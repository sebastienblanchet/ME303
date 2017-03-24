%Project #2: Part 3-4 - ME 303
%Sebastien Blanchet, Timothy Wulff

%Intialize script
close all
clear variables
clc

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

%Calculations
A = (pi/4)*(D_n^2);     % bar area [m^2]
S = P/(L*A);            % source strenght [W/m^3]
Q = S/(row*C_p);        % heat source term [deg.C/s]

%Define grid spacing
dx = 0.01;              % delta x [m]
dt = 1;                 % time step [s]

%Create difference and tolerance for grid independence
dif=1;
tol=0.1;
count=0;

%Create loop for iteration of grid independence
while dif>tol
    
    %Initialize the comparison matrix
    T_final=zeros(L/dx+1,t_end/dt+1);
    
    %Call the GetTnum function to perform CN for the first dt
    T_old_dt=GetTnum(dt);
    
    %Half the dt
    dt=dt-(dt/2);
    
    %Call the CN fxn to get the solution matrix for the half sized dt
    T_new_dt=GetTnum(dt);
    
    %Get every second value in the T_new_dt matrix and store in compare
    %matrix
    T_final(:,:)=T_new_dt(:,1:2:end);
    
    %Compare the different between temps for all respective times
    dif=max(max(abs(T_old_dt-T_final)));
    
    %Increase itteration count
    count=count+1;
    
    %As this is lengthy process, end if itterations exceed 5
    if count==5
        break
    end
    
end

%Call function to get the entire solution for the time spacing which converges
T_new_dt=GetTnum(dt);
x = 0:dx:L;        % x spacing vector
t = 0:dt:t_end;    % time vector

%3D Plot of temperature distribution
figure1=figure;
[x_3d,t_3d] = meshgrid(0:dx:L,0:dt:t_end);
surf(x_3d,t_3d,(T_new_dt)');
grid off
shading flat
xlabel('L [m]');
ylabel('t [s]');
zlabel('T [deg.C]');
title('Plot of t vs L vs T');

%Plot for desired times
plot_t=[1,(15/dt)+1,(30/dt)+1, (60/dt)+1, (120/dt)+1, (1000/dt)+1];
%Create table for all times
Table=T_new_dt(:,plot_t);
figure2=figure;
hold on
plot(x,Table);
xlabel('L [m]');
ylabel('T [deg.C]');
title('Plot of T vs L for various t');
legend('0s','15s','30s','60s','120s','1000s(SS)');
