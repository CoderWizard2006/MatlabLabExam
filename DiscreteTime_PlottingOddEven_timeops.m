
% DISCLAIMER--> ALL THESE CODES ARE COMPATIBLE WITH THE MATLAB 25 VERSION,
% not being compatible with matlab 2021a is not the author's problem, do
% not blame the author in case that happens.






function discrete_signal_analysis()

 
    n = -10:10;                          % discrete-time index
    x = @(n) sin(0.2*pi*n) + 0.5*cos(0.4*pi*n); % example discrete-time signal

  
    evenOdd(x, n);

    
    timeOperations(x, n);

end


function evenOdd(x, n)
    % Compute even and odd components
    x_even = (x(n) + x(-n))/2;
    x_odd  = (x(n) - x(-n))/2;

    % Plot
    figure;
    subplot(3,1,1);
    stem(n, x(n), 'b', 'filled'); hold on;
    title('Original Signal x[n]');
    xlabel('n'); ylabel('x[n]'); grid on;

    subplot(3,1,2);
    stem(n, x_even, 'r', 'filled');
    title('Even Component x_e[n]');
    xlabel('n'); ylabel('x_e[n]'); grid on;

    subplot(3,1,3);
    stem(n, x_odd, 'g', 'filled');
    title('Odd Component x_o[n]');
    xlabel('n'); ylabel('x_o[n]'); grid on;

end


function timeOperations(x, n)
    % Define transformations
    x_shift = @(n) x(n - 3);     
    x_rev   = @(n) x(-n);        
    x_scale = @(n) x(2*n);       

   
    figure;

    subplot(3,1,1);
    stem(n, x(n), 'b', 'filled'); hold on;
    stem(n, x_shift(n), 'r', 'filled');
    legend('x[n]','x[n-3]');
    title('Discrete-Time Shift Operation');
    xlabel('n'); ylabel('Amplitude'); grid on;

    subplot(3,1,2);
    stem(n, x(n), 'b', 'filled'); hold on;
    stem(n, x_scale(n), 'g', 'filled');
    legend('x[n]','x[2n]');
    title('Discrete-Time Scaling Operation');
    xlabel('n'); ylabel('Amplitude'); grid on;

    subplot(3,1,3);
    stem(n, x(n), 'b', 'filled'); hold on;
    stem(n, x_rev(n), 'm', 'filled');
    legend('x[n]','x[-n]');
    title('Discrete-Time Reversal Operation');
    xlabel('n'); ylabel('Amplitude'); grid on;

end





% Note that the only difference between the codes of discrete time and
% continous time signals is that the step is 1, and plot is replaced with
% stem.



function plotComplex()


    clc; clear; close all;


n = 0:20;  % Discrete-time index

% Example values for r and Omega
r = 1; 
Omega = 2; 

gamma = r/Omega;  % Real positive


% Case 1: 0 < gamma < 1 (decaying)

x1 = gamma.^n;


% Case 2: gamma > 1 (growing)

gamma2 = 2;  % gamma > 1
x2 = gamma2.^n;

% Case 3: complex gamma (oscillatory)

gamma3 = 0.9*exp(1j*pi/4);   % magnitude < 1, complex
x3 = gamma3.^n;
figure;

subplot(3,1,1);
stem(n, x1, 'filled');
title(['Decaying Exponential: \gamma = ' num2str(gamma)]);
xlabel('n'); ylabel('x[n]');
grid on;

subplot(3,1,2);
stem(n, x2, 'filled');
title(['Growing Exponential: \gamma = ' num2str(gamma2)]);
xlabel('n'); ylabel('x[n]');
grid on;

subplot(3,1,3);
stem(n, abs(x3), 'filled'); hold on;
stem(n, real(x3), 'r--'); 
stem(n, imag(x3), 'g--');
title(['Complex Exponential: \gamma = ' num2str(abs(gamma3)) 'e^{j\theta}']);
xlabel('n'); ylabel('x[n]');
legend('Magnitude','Real','Imag');
grid on;

sgtitle('Discrete-Time Exponential Signals x[n] = \gamma^n u[n]');



end