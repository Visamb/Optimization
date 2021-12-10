%% Armijo line search

function [lambda,ls_its] = armijo(f,alpha,epsilon)
% Line search method in accordance with Armijo's rule.
%
% :param f: input function
% :param alpha: double, alpha parameter for Armijo's rule.
% :param epsilon: double, epsilon parameter for Armijo's rule.
% :return: lambda: double, step length lambda.
% :return: ls_its: integer, number of iterations in the line search.

% Set the precision for the vpa
prec = 10;

%Ititialize count for linesearch iterations.
ls_its = 0;

%Increment for numerical approximation of derivatives.
h = vpa(1e-12,prec);

%Initialize lambda
lambda = vpa(1,64);

%Derivative at zer0 and test metrics according to Armijo's rule.
f_prime_0 = vpa((f(h)-f(0))/(h),prec);
T_lambda = f(0) + epsilon*lambda*f_prime_0;
T_alpha_lambda = f(0) + epsilon*alpha*lambda*f_prime_0;

%Terminate if the derivative in 0 s larger than 0.
if f_prime_0 > 0
    s = 'derivative larger than 0';
    return;
end

%If we are already in a minimum, terminate.
if f_prime_0 == 0
    lambda = 0;
    s = 'minima found in lambda=0'
end


while true
    
    %Increase counter
    ls_its = ls_its+1;
    
    %Test metrics according to Amrijo's rule.
    T_lambda = f(0) + epsilon*lambda*f_prime_0;
    T_alpha_lambda = f(0) + epsilon*alpha*lambda*f_prime_0;
  
    %Criterion: Termination criteria for nan values.
    if isnan(f(lambda))
        lambda = 0
        s = 'f(lambda) is not a number';
        break
    end
    
    %Criterion: If the derivative at lambda is zero, terminate.
    if abs(vpa((f(h+lambda)-f(lambda))/(h),prec)) < 1e-22
        s = 'really small derivative';
        break
    end
    
    %Criterion: For too small lambda, terminate loop.
    if lambda < 1e-22
        s = 'Lambda is really small';
        lambda = 0;
        break
    end
    
    %Counter for evaluation of Armijo's criteria.
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