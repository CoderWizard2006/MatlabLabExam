% DISCLAIMER--> ALL THESE CODES ARE COMPATIBLE WITH THE MATLAB 25 VERSION,
% not being compatible with matlab 2021a is not the author's problem, do
% not blame the author in case that happens.



function fourier()

clc; clear; close all;

To = pi;               % Fundamental period
wo = 2*pi/To;          % Fundamental angular frequency
h = 0.001;             % Time step
t = 0:h:(To-h);        % Time vector over one period

y = exp(-t/2);         % Example signal
N = length(y);

a0 = sum(y)/(N-1);


numTerms = 10;

a = zeros(1,numTerms);
b = zeros(1,numTerms);

for n = 1:numTerms
    a(n) = 2*sum(y.*cos(n*wo*t))/(N-1);
    b(n) = 2*sum(y.*sin(n*wo*t))/(N-1);
end

Cn = sqrt(a.^2 + b.^2);
thetan = atan2(-b,a); 


n = 0:numTerms;

figure;


subplot(2,2,1);
stem(n, [a0 a], 'k', 'filled');
xlabel('n'); ylabel('a_n');
title('Fourier Coefficients a_n');


subplot(2,2,2);
stem(n, [0 b], 'k', 'filled');
xlabel('n'); ylabel('b_n');
title('Fourier Coefficients b_n');

subplot(2,2,3);
stem(n, [a0 Cn], 'k', 'filled');
xlabel('n'); ylabel('C_n');
title('Amplitude Coefficients C_n');


subplot(2,2,4);
stem(n, [0 thetan], 'k', 'filled');
xlabel('n'); ylabel('\theta_n [rad]');
title('Phase Coefficients \theta_n');

sgtitle('Fourier Series Coefficients of y(t)');
grid on;


end


function manualDerived()

clc; clear; clf;
n= 1:10; an(1) = 0.504; an(n+1) = 0.504*2./(1+16*n.^2); %manually got an
bn(1) = 0; bn(n+1) = 0.504*8*n./(1+16*n.^2); %manually got bn
cn(1) = an(1);
cn(n+1) = sqrt(an(n+1).^2+bn(n+1).^2);
thetan(1) = 0; thetan(n+1) = atan2(-bn(n+1),an(n+1));
n=[0,n];
subplot(2,2,1); stem(n,an,'k'); ylabel('an'); xlabel('n');
subplot(2,2,2); stem(n,bn,'k'); ylabel('bn'); xlabel('n');
subplot(2,2,3); stem(n,cn,'k'); ylabel('cn'); xlabel('n');
subplot(2,2,4); stem(n,thetan,'k'); ylabel('\theta[rad]'); xlabel('n');

end




% To find the Dn
function exponentialFourier()

clc; clear; close all;


To = pi;                % Fundamental period
wo = 2*pi/To;           % Fundamental angular frequency
h = 0.001;              % Time step
t = 0:h:(To-h);         % One period


y = @(t) exp(-t/2);     % You can now do y(2*t), y(t-3), etc.

N = length(t);


numTerms = 10;               % max |n|
n = -numTerms:numTerms;      
Dn = zeros(size(n));        



% Note , if you want to plot the fourier coeff for f(at) or any other
% change, change y(t) to that.
for k = 1:length(n)
    Dn(k) = (1/To) * sum(y(t) .* exp(-1j*n(k)*wo*t)) * h; 
end


figure;

subplot(2,1,1);
stem(n, abs(Dn), 'k', 'filled');
xlabel('n'); ylabel('|D_n|');
title('Magnitude of Exponential Fourier Coefficients');
grid on;

subplot(2,1,2);
stem(n, angle(Dn), 'k', 'filled');
xlabel('n'); ylabel('Phase (rad)');
title('Phase of Exponential Fourier Coefficients');
grid on;

sgtitle('Exponential Fourier Series Coefficients D_n');


figure;
plot(t, y(t), 'b', 'LineWidth', 1.5); hold on;
plot(t, y(2*t), 'r--', 'LineWidth', 1.5); 
plot(t, y(t-1), 'g-.', 'LineWidth', 1.5);
xlabel('t'); ylabel('y(t)');
legend('y(t)','y(2t)','y(t-1)');
title('Original, Scaled, and Shifted Functions');
grid on;

end

exponentialFourier