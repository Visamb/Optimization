%% constrained_problem

%Starting point
x0 = [0.2;0.3;0.4;0.5;0.85];

%Chosen mu-parameters
mu_list = [1e1,1e2,1e3,1e4,1e5];

%Iterative solution. Choose specified problem and penalty function.
%The last returned optimal value from nonlinearmin is the found optimizer.
for i = mu_list
    mu = i;
    fprintf('%s\n', '------------------------------------------------------------------------------------------')
    fprintf('%s\n', 'Starting constrained optimization.')
    fprintf('%s %f\n', 'Current mu =',round(mu))
    fprintf('%s\n', '------------------------------------------------------------------------------------------')
    func = @(x) (sample_problem(x) + mu*h(x));
    [x0, no_its, normg] = nonlinearmin(func,x0,1e-6,1,0,1);
end

%% Functions

%Problem 9_5 (Problem C in project description)
function [y] = problem_9_5(x)

y = (x(1)-5)^2 + (x(2)-3)^2;

end

%Problem 9_3 (problem B in the problem description)
function [y] = problem_9_3(x)

y = exp(x(1)) + x(1)^2 + x(1)*x(2);

end

%Constrained problem A in the problem description
function [y] = sample_problem(x)

y = exp(x(1)*x(2)*x(3)*x(4)*x(5));

end

%Penalty function for constrained problem A
function [y] = h(x)

p = 2;

y = (x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 -10)^p + (x(2)*x(3)-5*x(4)*x(5))^p + (x(1)^3 + x(3)^3 +1)^p;

end

%Penalty function for constrained problem B
function [y] = h_9_3(x)

p = 2;

y = (0.5*x(1) + x(2) -1)^p;

end

%Penalty function for constrained problem C
function [y] = h_9_5(x)

p = 2;

% When g1,g2 < 0 we are in the right area.
% When g1,g2 > 0 we want the barrier to -> inf

g1 = max(0,(x(1)+x(2)-3));
g2 = max(0,(-x(1)+2*x(2)-4));

y = g1^p + g2^p;

end
