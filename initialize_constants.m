% Initialize constants for simulation

clear
clc
close all

%% Define COnstants

x0 = [85;
    0;
    0;
    0;
    0;
    0;
    0;
    0.1;
    0];


u = [0;
    -0.1;
    0;
    0.08;
    0.08];


TF = 60;


% RUNN

sim('RCAMsimulation.slx')


%% PLOTTING THE DATA

t = ans.simX.Time;



x1 = ans.simX.Data(:,1);
x2 = ans.simX.Data(:,2);
x3 = ans.simX.Data(:,3);
x4 = ans.simX.Data(:,4);
x5 = ans.simX.Data(:,5);
x6 = ans.simX.Data(:,6);
x7 = ans.simX.Data(:,7);
x8 = ans.simX.Data(:,8);
x9 = ans.simX.Data(:,9);






figure

subplot(3,3,1)
plot(t,x1)
legend('x_1')
grid on

subplot(3,3,4)
plot(t,x2)
legend('x_2')
grid on

subplot(3,3,7)
plot(t,x3)
legend('x_3')
grid on

subplot(3,3,2)
plot(t,x4)
legend('x_4')
grid on

subplot(3,3,5)
plot(t,x5)
legend('x_5')
grid on

subplot(3,3,8)
plot(t,x6)
legend('x_6')
grid on

subplot(3,3,3)
plot(t,x7)
legend('x_7')
grid on

subplot(3,3,6)
plot(t,x8)
legend('x_8')
grid on

subplot(3,3,9)
plot(t,x9)
legend('x_9')
grid on

