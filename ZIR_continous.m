% DISCLAIMER--> ALL THESE CODES ARE COMPATIBLE WITH THE MATLAB 25 VERSION,
% not being compatible with matlab 2021a is not the author's problem, do
% not blame the author in case that happens.



% question of finding ZIR might be given, find it using dsolve and then
% plot
function ZIR()

    syms y(t)
    Dy = diff(y, t);
    D2y = diff(y, t, 2);

  
    ode = D2y + 1.5*Dy + 0.5625*y == 0;

   
    cond1 = y(0) == 2;
    cond2 = Dy(0) == -3;

  
    sol(t) = dsolve(ode, cond1, cond2);

    disp("Zero Input Response Solution:");
    disp(sol(t));


    tn = 0:0.01:2*pi;

    sol_num = double(sol(tn));

 
    figure;
    plot(tn, sol_num, 'LineWidth', 1.5);
    title('Zero Input Response of the System');
    xlabel('t'); ylabel('y(t)');
    grid on;

end



% for plotting of h(t) , nothing will change except the inital conditions