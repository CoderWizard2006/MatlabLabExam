
function ZIR_discrete()
    clc; clear; close all;

    % -----------------------------
    % Example: Second-order difference equation
    % y[n+2] + 1.5*y[n+1] + 0.5625*y[n] = 0
    % with initial conditions y[0] = 2, y[1] = -3
    % -----------------------------

    % Coefficients of the difference equation
    a = [1 1.5 0.5625];   % characteristic equation coefficients
    b = [0 0 0];           % zero input, so RHS is 0

    % Initial conditions
    y0 = [2, -3];          % y[0] and y[1]

    % Number of points
    N = 20;
    y = zeros(1, N);

    % Set initial conditions
    y(1:length(y0)) = y0;

    % Compute ZIR recursively
    for n = length(y0)+1:N
        y(n) = -a(2)*y(n-1) - a(3)*y(n-2);  % difference eq rearranged
    end

    % Plot ZIR
    n_vec = 0:N-1;
    stem(n_vec, y, 'filled');
    title('Zero-Input Response (ZIR)');
    xlabel('n'); ylabel('y[n]');
    grid on;
end



function impulseResponse()

    clc; clear; close all;

% -----------------------------
% Discrete-time index
% -----------------------------
n = -2:10;  % time range
N = length(n);

% -----------------------------
% Initialize output y[n] with zero initial conditions
% -----------------------------
y = zeros(N,1);

% -----------------------------
% Define the input x[n] as a discrete impulse at n=0
% -----------------------------
x = zeros(N,1);
x(3) = 1;  % delta[n], since n(3) corresponds to n=0

% -----------------------------
% Iterative computation of y[n]
% Difference equation: y[n+2] - y[n+1] + 0.24*y[n] = x[n+2] - 2*x[n+1]
% -----------------------------
for k = 1:N-2
    y(k+2) = y(k+1) - 0.24*y(k) + x(k+2) - 2*x(k+1);
end

% -----------------------------
% Plot impulse response
% -----------------------------
stem(n, y, 'filled');
xlabel('n'); ylabel('h[n]');
title('Impulse Response h[n]');
grid on;



end





function ZSR()

    clc; clear; close all;


n = -2:10;  
N = length(n);

% -----------------------------
% Example system: difference equation
% y[n+2] - y[n+1] + 0.24*y[n] = x[n+2] - 2*x[n+1]
% Compute impulse response h[n]
% -----------------------------


y = zeros(N,1);
x_imp = zeros(N,1); 
x_imp(3) = 1; % delta[n] at n=0

for k = 1:N-2
    y(k+2) = y(k+1) - 0.24*y(k) + x_imp(k+2) - 2*x_imp(k+1);
end

h = y; % impulse response


x = [0 0 1 2 3 4 5]'; % e.g., n*u(n) or any discrete signal


y_out = conv(x, h); 
n_out = 0:length(y_out)-1; 

figure;

subplot(3,1,1);
stem(0:length(x)-1, x, 'filled');
title('Input Signal x[n]');
xlabel('n'); ylabel('x[n]');
grid on;

subplot(3,1,2);
stem(n, h, 'filled');
title('Impulse Response h[n]');
xlabel('n'); ylabel('h[n]');
grid on;

subplot(3,1,3);
stem(n_out, y_out, 'filled');
title('Output y[n] = x[n] * h[n]');
xlabel('n'); ylabel('y[n]');
grid on;




end



ZSR 