%Part C - ME 303 
%Sebastien Blanchet, Timothy Wulff

%Intialize code
close all
clear variables
clc

%Define C, lambda, intial delta r, tolerance and iteration counter
C = 4;
lambda = 3;
dr= 0.2;
tol = 0.0001; 
iterations=0;

%Declare symbolic variable x for plotting exact solution
syms x;

%Create void loop for grid spacing independence
while true
    
    %Start iteration counter
    iterations = iterations+1;
    
    %Generate  array for the r values from 0 to 10 with spapcing of dr
    r = 0:dr:10;
    r_len = length(r);

    %Initialize y and v arrays
    y = zeros(1,r_len);
    v = zeros(1,r_len);

    %Declare intial conditions given in Part B problem statement
    y(1)= C;
    v(1)= 0;
    
    %Create loop for Runge Kunta Calculations for all y and v array values
    for i= 1:r_len-1
        
        F_i = v(i);
        g_i = -v(i)-((lambda^2)*y(i)*r(i));
        
        %Account for 0/0 limit
        if g_i == 0 && r(i)== 0
            G_i = -18;
        else 
            G_i = g_i/r(i);
        end
        
        y_c =  y(i)+(dr/2)*F_i;
        v_c = v(i)+(dr/2)*G_i;
        F_c = v_c;
        G_c = -(v_c/(r(i)+(dr/2)))-((lambda^2)*y_c);
        y_cc = y(i)+(dr/2)*F_c;   
        v_cc = v(i)+(dr/2)*G_c;
        F_cc = v_cc; 
        G_cc =  -(v_cc/(r(i)+(dr/2)))-((lambda^2)*y_cc);
        y_s = y(i)+dr*F_cc; 
        v_s = v(i)+dr*G_cc;
        F_s = v_s;
        G_s = -(v_s/r(i+1))-((lambda^2)*y_s);
        y(i+1) = y(i)+dr*((F_i/6)+(F_c/3)+(F_cc/3)+(F_s/6));
        v(i+1) = v(i)+dr*((G_i/6)+(G_c/3)+(G_cc/3)+(G_s/6));
        
    end
       
    %Find exact values of y(r) for tolerance check
    y_exact = zeros(1,r_len);
    y_exact = 4*besselj(0,(3*r));   
    
    %Overwrite plotted data
    hold on
    
    %Plot y(r) wrt r from 0<=r<=10
    plot(r,y);
    
    %Check approximated and exact values of y with tolerances
    if abs(y-y_exact)<=tol
        %Get out of void loop if y values are converging
        %Plot exact bessel equation values for final dr values
        ezplot(4*besselj(0,(3*x)), [0, 10, -2, 4.1]);
        hold off
        break
    else
        %If y has not converged half the delta r and restart loop
        dr = dr-(dr/2);    
    end   
    
end

%Label x-axis as r
xlabel('r');

%Label y-axis as y
ylabel('y');

%Label plot title
title('Plot of y(r) wrt different dr');
