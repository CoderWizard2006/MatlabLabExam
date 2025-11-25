function reconstruct()


clc; clear; close all;

% -----------------------------
% Discrete-time signal parameters
% -----------------------------
N = 11;                 % number of points for one period
n = 0:N-1;              % discrete-time indices

% Example signal: x[n] = n*u[n] (one period)
x = n;

% -----------------------------
% Compute exponential Fourier coefficients Dn
% -----------------------------
Dn = zeros(1,N);
for k = 0:N-1
    Dn(k+1) = (1/N) * sum(x .* exp(-1j*2*pi*k*n/N));
end

% -----------------------------
% Reconstruct x[n] from Dn
% -----------------------------
x_recon = zeros(1,N);
for n_idx = 0:N-1
    x_recon(n_idx+1) = sum(Dn .* exp(1j*2*pi*(0:N-1)*n_idx/N));
end

% -----------------------------
% Plot original and reconstructed signals
% -----------------------------
figure;
stem(n, x, 'b', 'filled'); hold on;
stem(n, real(x_recon), 'r--', 'filled');
xlabel('n'); ylabel('Amplitude');
legend('Original x[n]', 'Reconstructed x[n]');
title('Reconstruction of x[n] from Exponential Fourier Series Coefficients');
grid on;



end







% DISCRETE FOURIER SERIES

function discreteFS()


clc; clear; close all;


%% Method 1: Direct DTFS Calculation

No = 32;                
n = 0:No-1;              

% Example signal: rectangular pulse + ones
xn = [ones(1,5) zeros(1,23) ones(1,4)];

% Preallocate DTFS coefficients
xr = zeros(1, No);

for r = 0:No-1
    xr(r+1) = sum(xn .* exp(-1j*r*2*pi/No*n)) / No;
end

% Plot magnitude and phase
figure('Name','DTFS: Direct Summation');
r = 0:No-1;
subplot(2,1,1);
stem(r, abs(xr), 'k', 'filled');
ylabel('|X[r]|'); xlabel('r');
title('Magnitude of DTFS Coefficients');

xr_round = round(1000*xr)/1000; % remove small computational errors
subplot(2,1,2);
stem(r, angle(xr_round), 'k', 'filled');
ylabel('Angle of X[r] (rad)'); xlabel('r');
title('Phase of DTFS Coefficients');



% USING FFT

xr_fft = fft(xn)/No;    % Normalize by period

figure('Name','DTFS: Using FFT');
subplot(2,1,1);
stem(r, abs(xr_fft), 'b', 'filled');
ylabel('|X[r]|'); xlabel('r');
title('Magnitude via FFT');

subplot(2,1,2);
stem(r, angle(xr_fft), 'b', 'filled');
ylabel('Angle of X[r] (rad)'); xlabel('r');
title('Phase via FFT');



end





% QUESTION TYPE WHERE DR coefficients are given


function discreteFS_example()

    clc; clear; close all;


% Given DTFS coefficients Dr

N = 12;
r = 0:N-1;

% Provided Fourier coefficients (Dr)
Dr = [0, -1j*0.8720, -1j*0.4330, 1j*0.3330, 1j*0.1440, -1j*0.2950,...
      0, 1j*0.2950, -1j*0.1440, -1j*0.3330, 1j*0.4330, 1j*0.8720];


% Reconstruct x[n] from Dr

n = -32:32;  % discrete-time indices (as given)
x = zeros(size(n));

for ni = 1:length(n)
    for ri = 1:N
        x(ni) = x(ni) + Dr(ri) * exp(1j*2*pi*(ri-1)*n(ni)/N);
    end
end


% Plot the reconstructed signal

figure;
stem(n, real(x), 'filled', 'b'); 
xlabel('n'); ylabel('x[n]');
title('Reconstructed Discrete-Time Signal from Given Fourier Coefficients');
grid on;



end








% finding the inverse fourier transform for a given signal


function inverseFTdiscrete()
    


    clc; clear; close all;

% Define frequency variable
Omega = linspace(-pi, pi, 1000);

% Define X(Omega) as given
X = zeros(size(Omega));
idx = (Omega >= -pi/2) & (Omega <= pi/2);
X(idx) = cos(Omega(idx));

% Compute inverse Fourier transform numerically
% Using definition:
% x(t) = (1/2pi) * integral X(Omega) * exp(j*Omega*t) dOmega
% over -pi to pi or wider frequency range

t = linspace(-10,10,1000);  % time vector

x = zeros(size(t));

for k = 1:length(t)
    integrand = X .* exp(1j * Omega * t(k));
    x(k) = (1/(2*pi)) * trapz(Omega, integrand);
end

% Plot real part of x(t)
figure;
plot(t, real(x), 'LineWidth', 2);
xlabel('t');
ylabel('x(t)');
title('Inverse Fourier Transform x(t)');
grid on;



end






% Fourier spectra changing with modifications


function simple_fourier_demo_loop()

    clc; clear; close all;

    % Define time vector
    t = linspace(-10,10,2000);
    dt = t(2) - t(1);

    % Original signal: rectangular pulse width 2
    x = @(t) double(abs(t) <= 1);

    % Frequencies to compute FT at
    omega = linspace(-20, 20, 500);

    % Preallocate FT results
    X = zeros(size(omega));
    Y1 = zeros(size(omega));
    Y2 = zeros(size(omega));
    Y3 = zeros(size(omega));

    % Define modified signals
    y1 = @(t) 3 * x(t);
    y2 = @(t) x(t - 2);
    y3 = @(t) x(t) .* exp(1j*5*t);

    % Compute FT for each omega using for-loop sum
    for k = 1:length(omega)
        w = omega(k);
        sum_x = 0; sum_y1 = 0; sum_y2 = 0; sum_y3 = 0;
        for n = 1:length(t)
            sum_x = sum_x + x(t(n)) * exp(-1j*w*t(n));
            sum_y1 = sum_y1 + y1(t(n)) * exp(-1j*w*t(n));
            sum_y2 = sum_y2 + y2(t(n)) * exp(-1j*w*t(n));
            sum_y3 = sum_y3 + y3(t(n)) * exp(-1j*w*t(n));
        end
        X(k) = sum_x * dt;
        Y1(k) = sum_y1 * dt;
        Y2(k) = sum_y2 * dt;
        Y3(k) = sum_y3 * dt;
    end

    % Plot time domain signals
    figure;
    subplot(3,1,1);
    plot(t, x(t), 'b', t, y1(t), 'r--', 'LineWidth', 1.5);
    title('Scalar Multiplication: x(t) and 3*x(t)');
    legend('x(t)', '3*x(t)');
    grid on;

    subplot(3,1,2);
    plot(t, x(t), 'b', t, y2(t), 'r--', 'LineWidth', 1.5);
    title('Time Shifting: x(t) and x(t-2)');
    legend('x(t)', 'x(t-2)');
    grid on;

    subplot(3,1,3);
    plot(t, x(t), 'b', t, real(y3(t)), 'r--', 'LineWidth', 1.5);
    title('Frequency Shifting: x(t) and x(t)*exp(j5t)');
    legend('x(t)', 'Modulated x(t)');
    grid on;

    % Plot magnitude spectra
    figure;
    subplot(3,1,1);
    plot(omega, abs(X), 'b', omega, abs(Y1), 'r--', 'LineWidth', 1.5);
    title('Magnitude Spectrum: Original and Scalar Multiplied');
    legend('X(\omega)', '3X(\omega)');
    grid on;

    subplot(3,1,2);
    plot(omega, abs(X), 'b', omega, abs(Y2), 'r--', 'LineWidth', 1.5);
    title('Magnitude Spectrum: Original and Time Shifted');
    legend('X(\omega)', 'e^{-j \omega 2}X(\omega)');
    grid on;

    subplot(3,1,3);
    plot(omega, abs(X), 'b', omega, abs(Y3), 'r--', 'LineWidth', 1.5);
    title('Magnitude Spectrum: Original and Frequency Shifted');
    legend('X(\omega)', 'Shifted Spectrum');
    grid on;

end


simple_fourier_demo_loop