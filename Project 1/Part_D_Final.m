%Part D - ME 303
%Sebastien Blanchet, Timothy Wulff

%Intialize code
close all
clear variables
clc

%Begining of code excerpt taken from Part C (lines 22 through 66)

%Define C and Lambda and intial spacing
C = 4;
lambda = 3;
dr= 0.1;
    
%Generate array for the r values from 0 to 10 with spapcing of dr
r = 0:dr:10;
r_len = length(r);

%Initialize y and v vectors
y = zeros(1,r_len);
v = zeros(1,r_len);

%Intial conditions
y(1)= C;
v(1)= 0;

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

%Exact values
y_exact = zeros(1,r_len);
y_exact = 4*besselj(0,(3*r));   

%end of Part C code excerpt

%Create empty arrays for numerical approximation and exact roots
r_num = [];
r_ex = [];

%Create symbolic variable to solve exact solution
syms X;

%Go through each value in y array
for i = 1 : r_len-1
    
    %Check if any value has crossed x-axis
    if y(i)/y(i+1) < 0
    
    %Interpolating for the numerical approximation
    r_num(end+1) = (-y(i)*dr/(y(i+1)-y(i)))+r(i); 
    
    %Find the exact value
    r_ex(end+1) = vpasolve(4*besselj(0,(3*X)) == 0 , X, [r(i), r(i+1)]);
    
    end

end

%Create an array for percentage error
root_num = length(r_num);
percentage = zeros(1,root_num);

%Go through each value in percentage array
for  j= 1 : root_num
    
    %Calculate difference between numerical and exact r's over exact x 100
    percentage(j)= ((r_num(j)-r_ex(j))/r_ex(j))*100;
    
end


    




