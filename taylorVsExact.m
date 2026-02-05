function taylorVsExact()
% taylorVsExact Compares Taylor series approximation with exact function values
%
%   taylorVsExact() interactively prompts the user for a mathematical function,
%   the number of Taylor series terms, and an approximation point. It then
%   generates two plots: one comparing the Taylor series approximation to the
%   exact function, and another showing the percentage error.
%
%   Interactive Inputs:
%       func - A mathematical function as a string (e.g., 'sin(x)', 'exp(-x)')
%              Special cases: use 'log(x)' for ln(x), 'log10(x)' for log base 10
%       n    - Number of terms for the Taylor series expansion (positive integer)
%       x0   - The center point for the Taylor series expansion (real number)
%
%   Outputs:
%       Two figures are generated:
%       1. 'TaylorVsExact' - Plots Taylor approximation (red dashed) vs exact
%                            function (black solid) over the interval [0, 6]
%       2. 'Error'         - Plots the percentage error between approximation
%                            and exact values
%
%   Examples:
%       taylorVsExact()
%       % Then enter at prompts:
%       %   Function: sin(x)
%       %   Number of terms: 5
%       %   Approximation point: 0
%
%   Mathematical Background:
%       The Taylor series expansion of f(x) around x0 is:
%       f(x) ≈ f(x0) + f'(x0)(x-x0) + f''(x0)(x-x0)^2/2! + ... + f^(n)(x0)(x-x0)^n/n!
%
%   Error Calculation:
%       MAE (Mean Absolute Error) = |exact - approximation|
%       TRE (True Relative Error) = (exact - approximation) / exact * 100%
%
%   Requirements:
%       - Symbolic Math Toolbox
%       - Function must be differentiable n times at x0
%
%   See also: TAYLOR, DIFF, SYM, SUBS
%
%   Author: Rebanta-Gupta
%   Repository: https://github.com/Rebanta-Gupta/taylor-vs-error

% -------------------------------------------------------------------------

func = input("Please enter an acceptable function (ex. sin(x), exp(-x)), also for ln(x) use log(x) and for normal log(x) use log10(x): ", "s");
n = input("Please enter the Number of Terms for the Taylor Series: ");
x0 = input("Please enter Approximate Point, x0: ");

x = sym("x");

x_plot = linspace(0, 6, 50);

taylor_pts = sym(zeros(1, (n+1)));

try
    f = str2sym(func);
    f_x0 = double(subs(f, x, x0));
    taylor_pts(1) = f_x0;
catch
    error("%s is not an acceptable function!", func);
end

y_estimate_vals = my_taylor(x_plot, x0, n, f, taylor_pts);
y_exact = subs(f, x, x_plot);

error_mae = y_exact - y_estimate_vals;
MAE = double(mean(abs(error_mae(~isnan(error_mae))), "all"));
error_tre = error_mae ./ y_exact;
abs_error = abs(error_tre);
error_tay = abs_error * 100;

for a = 1:50
    if error_tay(a) < 1e-6
        error_tay(a) = 0;
    end
end

figure("Name", "TaylorVsExact", "NumberTitle","off")
hold on 
title(func+" from taylor series expression compared to exact centered around "+num2str(x0))
xlabel('x')
ylabel(func)
plot(x_plot, y_estimate_vals, 'r--','Displayname','Taylor');
plot(x_plot, y_exact, 'k-','Displayname','Exact');
grid on
legend show

figure("Name", "Error Analysis", "NumberTitle","off")
hold on 
title("Error "+func+" from taylor series expression centered around "+num2str(x0))
xlabel('x')
ylabel("True Relative Error")
plot(x_plot,error_tay,'b-.');
annotation('textbox', [0.75, 0.75, 0.3, 0.15], ...
           'String', "MAE = "+num2str(MAE), ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', ...
           'EdgeColor', 'black', ...
           'FontSize', 10, ...
           'FontWeight', 'bold');
hold off

    function y_vals = my_taylor(x_plot, x0, n, f, taylor_pts)
        y_vals = zeros(1, 50);
        for k = 1:50
            h = (x_plot(k) - x0);
            df = sym(zeros(1,n));
            for i = 1:n
                df(i) = diff(f, i);
            end
            y_df = double(subs(df, x, x0));
            fact = factorial(1:n);
            numerator = sym(zeros(1, n));
            for i = 1:n
                numerator(i) = (h^(i)) * y_df(i);
                taylor_pts(i+1) = numerator(i) / fact(i);
            end
            y_vals(k) = sum(taylor_pts, "all");
        end
    end
end

% -------------------------------------------------------------------------
