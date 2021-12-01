%% constrained_problem

%Starting point
x0 = [0;0];

%Chosen mu-parameters
mu_list = vpa([1, 1e1, 1e2, 1e3, 1e4, 1e7, 1e8, 1e9],64);

%Iterative solution. Choose specified problem and penalty function.
for i = mu_list
    mu = i;
    fprintf('%s\n', '------------------------------------------------------------------------------------------')
    fprintf('%s\n', 'Starting constrained optimization.')
    fprintf('%s %f\n', 'Current mu =',round(mu))
    fprintf('%s\n', '------------------------------------------------------------------------------------------')
    func = @(x) (problem_9_5(x) + mu*h_9_5(x));
    [x0, no_its, normg] = nonlinearmin(func,x0,1e-6,1,0,1);
end

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


function [y] = h_9_5(x)

p = 2;

% When g1,g2 < 0 we are in the right area.
% When g1,g2 > 0 we want the barrier to -> inf

g1 = max(0,(x(1)+x(2)-3));
g2 = max(0,(-x(1)+2*x(2)-4));

y = g1^p + g2^p;

end
