clc; clear;

Age = 22;
x0 = 0;
xfinal = 1 + 0.01 * Age;
N_values = [2, 5, 10, 20, 40];

exactSol1 = @(x) 7*exp(x) - x.^3 - 3*x.^2 - 6*x - 6;
exactSol2 = @(x) (70/9)*exp(-0.3*x) - (43/9)*exp(-1.2*x);

f1 = @(x, y) x.^3 + y;
f2 = @(x, y) -1.2*y + 7*exp(-0.3*x);

methods = {'RK2', 'RK3', 'RK4'};

errors1 = zeros(length(N_values), 3); 
errors2 = zeros(length(N_values), 3);

fprintf('Solutions for each step:\n');

for idx = 1:length(N_values)
    N = N_values(idx);
    h = (xfinal - x0) / N;
    x = x0:h:xfinal;

    fprintf('\n>> N = %d, h = %.4f\n', N, h);
    
    for m = 1:3
        method = methods{m};
        
        fprintf('\nMethod: %s - Question 1\n', method);
        y1 = solve_ODE(f1, x0, 1, h, N, method);
        y1_exact = exactSol1(x(end));
        errors1(idx, m) = abs(y1(end) - y1_exact);
        fprintf('Final y value: %.8f - Exact: %.8f - Error: %.2e\n', y1(end), y1_exact, errors1(idx, m));
        
        fprintf('\nMethod: %s - Question 2\n', method);
        y2 = solve_ODE(f2, x0, 3, h, N, method);
        y2_exact = exactSol2(x(end));
        errors2(idx, m) = abs(y2(end) - y2_exact);
        fprintf('Final y value: %.8f - Exact: %.8f - Error: %.2e\n', y2(end), y2_exact, errors2(idx, m));
    end
end

figure;
for m = 1:3
    loglog((xfinal-x0)./N_values, errors1(:, m), '-o', 'DisplayName', methods{m}); hold on;
end
title('Q1 Error vs step size');
xlabel('Step size h'); ylabel('Absolute error');
legend('Location', 'southeast'); grid on;

figure;
for m = 1:3
    loglog((xfinal-x0)./N_values, errors2(:, m), '-o', 'DisplayName', methods{m}); hold on;
end
title('Q2 Error vs step size');
xlabel('Step size h'); ylabel('Absolute error');
legend('Location', 'southeast'); grid on;

fprintf('\n here are the comments \n');
fprintf('RK4 almost always had the lowest error rate, with RK3 and RK2 following.\n');
fprintf('RK4 rapidly improved, had the lowest error and was most efficient overall.\n');

function y = solve_ODE(f, x0, y0, h, N, method)
    y = zeros(1, N+1);
    y(1) = y0;
    x = x0;

    for i = 1:N
        switch method
            case 'RK2' 
                k1 = f(x, y(i));
                k2 = f(x + h, y(i) + h * k1);
                y(i+1) = y(i) + (h/2)*(k1 + k2);
            case 'RK3'
                k1 = f(x, y(i));
                k2 = f(x + h/2, y(i) + h/2 * k1);
                k3 = f(x + h, y(i) - h * k1 + 2 * h * k2);
                y(i+1) = y(i) + (h/6)*(k1 + 4*k2 + k3);
            case 'RK4'
                k1 = f(x, y(i));
                k2 = f(x + h/2, y(i) + h/2 * k1);
                k3 = f(x + h/2, y(i) + h/2 * k2);
                k4 = f(x + h, y(i) + h * k3);
                y(i+1) = y(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
        x = x + h;
        fprintf('Step %2d: x = %.4f, y = %.8f\n', i, x, y(i+1));
    end
end
