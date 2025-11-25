function oddEven()

clear; clc; close all;

t = linspace(-5,5,2000);     


x = @(t) sin(2*t) + 0.5*cos(3*t) + exp(-abs(t));


x_t  = x(t);
x_nt = x(-t);


x_even = (x_t + x_nt)/2;
x_odd  = (x_t - x_nt)/2;


figure;

subplot(3,1,1);
plot(t, x_t, 'LineWidth', 1.5);
title('Original Continuous-Time Signal x(t)');
xlabel('t'); ylabel('x(t)'); grid on;

subplot(3,1,2);
plot(t, x_even, 'LineWidth', 1.5);
title('Even Component x_e(t)');
xlabel('t'); ylabel('x_e(t)'); grid on;

subplot(3,1,3);
plot(t, x_odd, 'LineWidth', 1.5);
title('Odd Component x_o(t)');
xlabel('t'); ylabel('x_o(t)'); grid on;

end



function timeshifting() 
    


    t=0:0.01:2*pi;

    f= @(t) sin(t) + cos(t);

    f_t = f(t);
    plot(t, f_t, 'LineWidth', 1.5);
    title('Original Function f(t)');
    xlabel('t'); ylabel('f(t)'); grid on;hold on;

    f2_t= f(2*t);
    plot(t,f2_t,'LineWidth',1.5);
    title('original function f(t)');
    
   % similarly all other operations of shifting can be done.

end





% sketching of functions which are periodic 
function trapezoid_clean()

clear; clc; close all;


t = linspace(-6, 6, 2000);


x = max(0, min(1, 5 - abs(t)));


figure;
plot(t, x, 'LineWidth', 1.5);
title('Trapezoidal Pulse x(t)');
xlabel('t'); ylabel('x(t)'); grid on;


energy = trapz(t, x.^2);
fprintf("Total Energy = %.4f\n", energy);


T = 10;
t2 = linspace(-15, 15, 4000);


tp = mod(t2 + 5, T) - 5;

x_periodic = max(0, min(1, 5 - abs(tp)));


figure;
plot(t2, x_periodic, 'LineWidth', 1.5);
title('Periodic Trapezoid of Period 10 (3 cycles)');
xlabel('t'); ylabel('x_p(t)'); grid on;

end


trapezoid_clean