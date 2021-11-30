%% constrained_problem

x0 = [-2;2;2;-1;-1];
mu_list = vpa([1e1, 1e2, 1e3, 1e4, 1e5, 1e6],64)

for i = mu_list
    current_mu = i;
    func = @(x) (sample_problem(x) + i*h(x));
    [x0, no_its, normg] = nonlinearmin(func,x0,1e-3,1,0);
end

minimum = x0;



%%
test = h(x0)

%% Functions
function [y] = sample_problem(x)

y = exp(x(1)*x(2)*x(3)*x(4)*x(5));

end

function [y] = h(x)

p = 4;

y = (x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 -10)^p + (x(2)*x(3)-5*x(4)*x(5))^p + (x(1)^3 + x(3)^3 +1)^p;

end

