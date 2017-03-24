%intialize code
close all
clear variables
clc

%Generate an array for the t values from 0 to 2 with spapcing of 0.01
deltat=0.01;
t= 0:deltat:2;
colnum=size(t, 2);

%Call f(t) to find all values of fxn for delta t
f_t = ((4/5)*cos(2*t))+((11/10)*sin(2*t))+((1/5)*exp(-t));

%Creates and an array the same length as f_t consisting of only zeros to
%speed up code
df_dt=zeros(1,colnum);

%go through each
for i= 1:colnum
    
    %Foward differencing at t = 0
    if t(i) == 0
        df_dt(i)=(f_t(i+1)-f_t(i))/(deltat);
    
    %Backward differencing at t = 2
    elseif t(i) == 2
        df_dt(i)=(f_t(i)-f_t(i-1))/(deltat);
    
    %Central differencing everywhere else
    else
        df_dt(i)=(f_t(i+1)-f_t(i-1))/(2*deltat);
    end

end

%Plot f(t) and df/dt wrt t from 0<=t<=2
plot(t,f_t, t,df_dt);

%Label x axis as t
xlabel('t');

%Label y axis as y
ylabel('y');

%Label plot title
title('Plot of f(t) and df/dt');

%Create legend to differentiate between data sets
legend('f(t)','df/dt');

%Select arrays to export
CSV_file=[t', f_t', df_dt'];
    
%Export the .CSV file
csvwrite('Part_A_Export.csv',CSV_file);

        