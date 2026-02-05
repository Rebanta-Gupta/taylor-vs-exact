function taylorVsExact()
% This function plots a graph of the provided function with the 
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
error_tre = error_mae ./ y_exact;
abs_error = abs(error_tre);
error_tay = abs_error * 100;

figure("Name", "TaylorVsExact", "NumberTitle","off")
hold on 
title(func+" from taylor series expression compared to exact centered around "+num2str(x0))
xlabel('x')
ylabel(func)
plot(x_plot, y_estimate_vals, 'r--','Displayname','taylor');
plot(x_plot, y_exact, 'k-','Displayname','exact');
grid on
legend show

figure("Name", "Error", "NumberTitle","off")
hold on 
title("error "+func+" from taylor series expression centered around "+num2str(x0))
xlabel('x')
ylabel("error")
plot(x_plot,error_tay,'b-.');
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