%% linesearch

function [lambda,nbr_it,history] = linesearch(f,tol)

%Initialize number of iterations
nbr_it = 0;

%Max number of iterations
max_it = 100000;

%Numerical precision for floats
prec = 32;

%Initialize a and b for bisection method
a = 0;
b = vpa(1e2,prec);

epsilon = vpa(1e-18,prec);

history = [];
       
while b-a > tol && nbr_it < max_it
    
    %Find point lambda_k and evaluate the derivative in lambda_k
    lambda = (a+b)/2;
    f_prime = vpa((f(lambda+epsilon)-f(lambda))/(epsilon),prec);
        
    %Check condition and set new a,b.
    if f_prime <= 0
        a = lambda;
    else
        b = lambda;
    end
    
    %Update nbr and history
    nbr_it = nbr_it + 1;
    history = [history; a,b, b-a];
end
    
    
lambda = lambda;
    
    
    
    
    

    
   % history = [history f(lambda)];
    
   % if abs(f_prime) < tol || abs(vpa(f(lambda),prec)-old_f)/abs(old_f) < tol
    %    abs(vpa(f(lambda),prec)-old_f)/abs(old_f)
     %   nbr_it = i;
     %   break
    %end
%end

if isnan(f(lambda)) || f(lambda)>f(0)
error('Bad job of the line search!')
end


end