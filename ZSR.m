
% DISCLAIMER--> ALL THESE CODES ARE COMPATIBLE WITH THE MATLAB 25 VERSION,
% not being compatible with matlab 2021a is not the author's problem, do
% not blame the author in case that happens.



function question1a()

    % Define signals
    x = @(t) (1./(t.^2 +1));                       % Unit step
    h = @(t) (1 ) .* (t >= 0);      % Example impulse response

    % Time vectors
    tvec = -0.25:0.1:5;     % this can be large too
    dtau = 0.1;            
    tau = -3.5:dtau:5;       

    
    y = zeros(1,length(tvec));

    % Counter
    ti = 0;

    % Convolution loop with animation
    for t = tvec
        ti = ti + 1;

        % Compute integrand x(τ) * h(t-τ)
        xh = x(tau) .* h(t - tau);

        % Numerical integration
        y(ti) = trapz(tau, xh);

        % Plot signals
        figure(1); clf;
        subplot(2,1,1);
        plot(tau, x(tau), 'r', 'LineWidth', 1.5); hold on;
        plot(tau, h(t - tau), 'm', 'LineWidth', 1.5);
        plot(t, 0, 'ro', 'MarkerFaceColor','r');
        xlabel('\tau'); ylabel('Amplitude');
        legend('x(\tau)','h(t-\tau)','Location','Best');
        title(['Overlap at t = ', num2str(t)]);

        subplot(2,1,2);
        plot(tvec, y, 'b', 'LineWidth', 1.5); hold on;
        plot(tvec(ti), y(ti), 'ro', 'MarkerFaceColor','r');
        xlabel('t'); ylabel('y(t)');
        title('Convolution Output');
        drawnow; pause(0.01);
    end

    % Final plot (clean)
    figure(2);
    plot(tvec, y, 'b', 'LineWidth', 2);
    xlabel('t'); ylabel('y(t)');
    title('Final Convolution Result');
    grid on;


end

%NOTE->

% only the function given in question will be changed, and take the largest
% possible time period as vector, even if it is extra that convolution will
% be zero.





% OR alternatively using dsolve and conv function
% finding ZSR if a equation is given

function ZSR_with_convolution_corrected()

    
    syms y(t) t
    Dy = diff(y, t);
    D2y = diff(y, t, 2);

   
    ode = D2y + 1.5*Dy + 0.5625*y == dirac(t);
    conds = [y(0) == 0, Dy(0) == 1];  

    % Solve for impulse response h(t)
    hSol(t) = dsolve(ode, conds);

    disp('Impulse Response h(t):');
    pretty(hSol);


    t_num = linspace(0, 10, 1000);
    h_num = double(hSol(t_num));

    x = @(t) exp(-t) .* (t >= 0);

    x_num = x(t_num);

    % Convolution
    dt = t_num(2) - t_num(1);
    y_num = conv(x_num, h_num) * dt;

    % Time vector for convolution output
    t_conv = linspace(t_num(1)+t_num(1), t_num(end)+t_num(end), length(y_num));

    % Plot results
    figure;
    subplot(3,1,1);
    plot(t_num, x_num, 'LineWidth', 1.5); title('Input x(t)'); grid on;
    subplot(3,1,2);
    plot(t_num, h_num, 'LineWidth', 1.5); title('Impulse Response h(t)'); grid on;
    subplot(3,1,3);
    plot(t_conv, y_num, 'LineWidth', 1.5); title('Output y(t) = x(t) * h(t)'); grid on;

end








% below is the matlab code for convolution of discrete time signals ( it is
% the modified version of the ZSR of continous time signals

function ZSR_discrete()

    n = 0:20;                  % discrete-time index
    x = @(n) 1./(n.^2 + 1);    % example input signal
    h = @(n) 1 .* (n >= 0);    % example impulse response

    y = zeros(1,length(n));   


    for ni = 1:length(n)
        k = 0:n(ni);           % summation index from 0 to current n
        xh = x(k) .* h(n(ni)-k); % compute product
        y(ni) = sum(xh);       % discrete summation

       
        figure(1); clf;

        subplot(2,1,1);
        stem(k, x(k), 'r', 'filled'); hold on;
        stem(k, h(n(ni)-k), 'm', 'filled');
        xlabel('k'); ylabel('Amplitude');
        legend('x[k]','h[n-k]','Location','Best');
        title(['Overlap at n = ', num2str(n(ni))]);

        subplot(2,1,2);
        stem(n, y, 'b', 'filled'); hold on;
        stem(n(ni), y(ni), 'ro', 'filled');
        xlabel('n'); ylabel('y[n]');
        title('Convolution Output');
        drawnow; pause(0.05);
    end


    figure(2);
    stem(n, y, 'b', 'filled');
    xlabel('n'); ylabel('y[n]');
    title('Final Convolution Result (Discrete-Time)');
    grid on;

end


ZSR_discrete