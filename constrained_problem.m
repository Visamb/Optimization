%% constrained_problem

x0 = [0;0];
mu_list = vpa([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7],64)

for i = mu_list
    current_mu = i;
    func = @(x) (problem_9_5(x) + i*g_9_5(x));
    [x0, no_its, normg] = nonlinearmin(func,x0,1e-5,1,0);
end

minimum = x0;



%%
test = h(x0)

%% Functions


function [y] = problem_9_5(x)

y = (x(1)-5)^2 + (x(2)-3)^2;

end

function [y] = problem_9_3(x)

y = exp(x(1)) + x(1)^2 + x(1)*x(2);

end


function [y] = sample_problem(x)

y = exp(x(1)*x(2)*x(3)*x(4)*x(5));

end

function [y] = h(x)

p = 2;

y = (x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 -10)^p + (x(2)*x(3)-5*x(4)*x(5))^p + (x(1)^3 + x(3)^3 +1)^p;

end

function [y] = h_9_3(x)

p = 2;

y = (0.5*x(1) + x(2) -1)^p;

end

function [y] = g_9_5(x)

p = 1;

y = (3-x(1)-x(2))^-p + (4+x(1)-2*x(2))^-p;

end
