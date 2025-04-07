clear all; clc; close all;

%% Task set A 
%% 2.1

%% 2.1.2
P = [0.7, 0.4, 0,   0.2;   % Row 1: Next state = S
     0.3, 0,   0,   0;     % Row 2: Next state = E
     0,   0.3, 0,   0;     % Row 3: Next state = I
     0,   0.3, 1.0, 0.8];  % Row 4: Next state = R

%% 2.1.3
% Present state vector
x = [0, 1, 0 ,0]';

% Next state vector
x_next = P * x;

% Display the result
disp('2.1.3');
disp(x_next);

%% 2.1.4
clear x, x_next;

% Present state vector
x = [1, 0, 0, 0]';

% 5 steps state vector
x_next = P^5 * x;

% Display the result
disp('2.1.4');
disp(x_next);

%% Task set B
clear x, x_next;
%% 3.1.1.a

x = zeros(4, 32); 
x(:,1) = [1; 0; 0; 0]; 

for n = 2:32
    x(:,n) = P * x(:,n-1);
    x(:,n) = x(:,n) / norm(x(:,n), 1);  % normalize
end

% Graph making
figure (1);
days = 0:31;

plot(days, x(1,:), '-o', days, x(2,:), '-x', days, x(3,:), '-s', days, x(4,:), '-^');
legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
xlabel('Day');
ylabel('Probability');
title('SEIR Model - State Probabilities Over Time');

%% 3.1.2.a

y = zeros(4, 32); 
y(:,1) = [0.15; 0.85; 0; 0]; 

for n = 2:32
    y(:,n) = P * y(:,n-1);
    y(:,n) = y(:,n) / norm(y(:,n), 1);  % normalize
end

% Graph making
figure (2);

plot(days, y(1,:), '-o', days, y(2,:), '-x', days, y(3,:), '-s', days, y(4,:), '-^');
legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
xlabel('Day');
ylabel('Probability');
title('SEIR Model - State Probabilities Over Time');

%% 3.1.3
% Calculate the difference in the first 10 steps, then graph

for n = 1:10
    diff(:,n) = x(:,n) - y(:,n);
end

figure (3);

plot(1:10, diff(1,:), '-o', 1:10, diff(2,:), '-x', 1:10, diff(3,:), '-s', 1:10, diff(4,:), '-^');
legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
xlabel('Day');
ylabel('Difference');
title('Difference in State Probabilities Over Time');

%% 3.4.a

[V, D] = eig(P);

% Display the result
disp('4.a');
disp(V);

%% 3.4.d

lambda = diag(D);
index = 3;
x_stationary = V(:,index);
x_stationary = x_stationary / norm(x_stationary, 1);

% Display the result
errors = zeros(1, 32);

for n = 1:32
    errors(n) = norm(x_stationary - x(:,n), 1);
end

% Graph making
figure (4);
semilogy(0:31, errors, 'b-o', 'LineWidth', 1.5);
xlabel('Iteration Set (n)');
ylabel('Absolute Error');
title('Convergence to stationary distribution');

% Shows us that the convergence is exponential (linear, but on a log scale)

%% 3.5.b

Pnew = [0.7, 0.4, 0,   0.1, 0;   % Row 1: Next state = S
        0.3, 0,   0,   0,   0;     % Row 2: Next state = E
        0,   0.3, 0,   0,   0;     % Row 3: Next state = I
        0,   0.3, 1.0, 0.8, 0;  % Row 4: Next state = R
        0,   0,   0,   0.1,   1];  % Row 5: Next state = Im

%% 3.5.c

x = zeros(5, 251); 
x(:,1) = [1; 0; 0; 0; 0]; 

for n = 2:251
    x(:,n) = Pnew * x(:,n-1);
    x(:,n) = x(:,n) / norm(x(:,n), 1);  % normalize
end
        
% Graph making
figure (5);
days = 0:250;

plot(days, x(1,:), '-', days, x(2,:), '-', days, x(3,:), '-', days, x(4,:), '-', days, x(5,:), '-', 'LineWidth', 1.5);
legend('Susceptible', 'Exposed', 'Infected', 'Recovered', 'Immunized');
xlabel('Day');
ylabel('Probability');
title('SEIR Model - State Probabilities Over Time');

% Well ofc its different... its a new model!

%% Task Set C
clear all; 
%% 4.1.1

P = [
    0.7, 0,   0,  0.2, 0,    0;   % S
    0.3, 0,   0,  0,   0,    0;   % E
    0,   0.5, 0,  0,   0,    0;   % I
    0,   0.5, 1,  0.8, 0,    0;   % R
    0,   0,   0,  0,   0.25, 0;   % V
    0,   0,   0,  0,   0.75, 1.0  % Im
];


%% 4.1.2

disp('4.1.2');
[V, D] = eig(P);
eigenvalues = real(D);  % Get only the real parts

disp(eigenvalues);
% Multiplicty is 2 for lambda = 1

%% 4.1.3

x = zeros(6, 31);
x(:,1) = [1; 0; 0; 0; 0; 0]';

for n = 2:31
    x(:,n) = P * x(:,n-1);
    x(:,n) = x(:,n) / norm(x(:,n), 1);  % normalize
end

% Graph making
figure (6);
days = 0:30;

plot(days, x(1,:), '-o', days, x(2,:), '-x', days, x(3,:), '-s', days, x(4,:), '-^', days, x(5,:), '-', days, x(6,:), '-', 'LineWidth', 1.5);
legend('Susceptible', 'Exposed', 'Infected', 'Recovered', 'Vaccinated', 'Immunized');
xlabel('Day');
ylabel('Probability');
title('SEIR-VIm Model - State Probabilities Over Time');

% Exact same distribution as if there was no vaccination / immune state

%% 4.1.4

x = zeros(6, 31);
x(:,1) = [.33; 0; 0; 0; .67; 0]';

for n = 2:31
    x(:,n) = P * x(:,n-1);
    x(:,n) = x(:,n) / norm(x(:,n), 1);  % normalize
end

% Graphe making
figure (7);

plot(days, x(1,:), '-o', days, x(2,:), '-x', days, x(3,:), '-s', days, x(4,:), '-^', days, x(5,:), '-', days, x(6,:), '-', 'LineWidth', 1.5);
legend('Susceptible', 'Exposed', 'Infected', 'Recovered', 'Vaccinated', 'Immunized');
xlabel('Day');
ylabel('Probability');

%% Saving the Figures
% numFigures = 7;
% for i = 1:numFigures
%     fig = figure(i);
%     % pause(5)
%     print(fig, sprintf('/Users/rowdyer/Documents/Coding/APPM2360_Project_2/Figures/figure%d', i), '-dpng', '-r1000');
% end