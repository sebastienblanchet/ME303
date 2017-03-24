%This function takes in a time step, performs Crank Nicolson and returns
%the numerical solution for the temperature distribution

function [Tn] = GetTnum(dt)

%Define constants
L = 0.15;               % length [m]
D_n = 0.01125;          % diameter [m]
alpha = 1.17e-4;        % thermal diffusivity [m^2/s]
P = 8;                  % input power [W]
row = 8933;             % density [kg/m^3]
C_p = 385;              % specific heat capacity [J/kg*deg.C]
K = 401;               % thermal conductivity [W/m]
T_0 = 18;               % outside temp. [deg.C]
T_L = 25;               % final temp. [deg.C]
t_end = 1000;           % end time [s]
dx = 0.01;              % delta x [m]

%Calculations
A = (pi/4)*(D_n^2);     % bar area [m^2]
S = P/(L*A);            % source strenght [W/m^3]
Q = S/(row*C_p);        % heat source term [deg.C/s]
F = (alpha*dt)/(dx^2);  % Fouriers number []

%Create vectors and matrices
To= zeros(L/dx+1,t_end/dt+1);   % old temperature

%Define IC & Dirichlet BC
To(:,1) = T_0;                  %T(x,0)= T0 for all x
To(L/dx+1,:) = T_L;             %T(L,t)= TL for all t

%Initialize solution matrix for new temperature
Tn=To;

%Create tolerance for itteration convergence
tolerance=0.00001;
convergence=1;

%Create itteration loop for convergence of temperature matrix
while (tolerance<convergence)
    for k=1:t_end/dt
        %Discretized Neumann BC
        Tn(1,k+1)=Tn(2,k+1);    %dT/dx(0,t)=0 for all t
        for i=2:L/dx
            %Crank-Nicolson equation
            Tn(i,k+1) = ((((1-F)*To(i,k))+((F/2)*(To(i+1,k)+To(i-1,k)+Tn(i+1,k+1)+Tn(i-1,k+1)))+(Q*dt))/(1+F));
        end
    end
    
    %Max difference between new and old temp., set equal to convergence.
    convergence=max(max(abs(Tn-To)));
    To=Tn;
  
end

end

