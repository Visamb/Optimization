%% Armijo line search

function [lambda,ls_its] = armijo2(f,alpha,epsilon,f_prime_0)
% Line search method in accordance with Armijo's rule.
%
% :param f: input function
% :param alpha: double, alpha parameter for Armijo's rule.
% :param epsilon: double, epsilon parameter for Armijo's rule.
% :return: lambda: double, step length lambda.
% :return: ls_its: integer, number of iterations in the line search.

%Ititialize count for linesearch iterations.
ls_its = 0;

%Increment for numerical approximation of derivatives.
h = 1e-12;

%Initialize lambda
lambda = 1;

%Derivative at zer0 and test metrics according to Armijo's rule.
%f_prime_0 = (f(h)-f(0))/(h);

%Terminate if the derivative in 0 s larger than 0.
if f_prime_0 > 0
    s = 'derivative larger than 0';
    return;
end

%If we are already in a minimum, terminate.
if f_prime_0 == 0
    lambda = 0;
    s = 'minima found in lambda=0';
end

while true
    
    %Increase counter
    ls_its = ls_its+1;
    
    %Test metrics according to Amrijo's rule.
    f_0 = f(0);
    T_lambda = f_0 + epsilon*lambda*f_prime_0;
    T_alpha_lambda = f_0 + epsilon*alpha*lambda*f_prime_0;
    
    %Criterion: For too small lambda, terminate loop.
    if lambda < 1e-22
        lambda = 0;
        break
    end
    
    num = 0;
    
    %Armijo's criteria.
    if f(lambda) > T_lambda
        lambda = lambda/alpha;
        num = num+1;
    end
    
    if f(lambda*alpha) < T_alpha_lambda 
        lambda = lambda*alpha;
        num = num+1;
    end
        
    %If none of the inequalities were fulfilled, break the loop and end
    %search.
    if num == 0
        break
    end
end

if isnan(f(lambda)) || f(lambda)>f(0)
error('Bad job of the line search!')
end

end