%% armijo line search

function [lambda,ls_its] = armijo(f,alpha,epsilon)

prec = 90;
ls_its = 0;
h = vpa(1e-12,prec);

lambda = vpa(3,64);
f_prime_0 = vpa((f(h)-f(0))/(h),prec);
T_lambda = f(0) + epsilon*lambda*f_prime_0;
T_alpha_lambda = f(0) + epsilon*alpha*lambda*f_prime_0;

if f_prime_0 > 0
    s = 'derivative larger than 0';
    return;
end

if f_prime_0 == 0
    lambda = 0;
    s = 'minima found in lambda=0'
end




%;
%vpa((f(h+lambda)-f(lambda))/(h),prec)
while true
    ls_its = ls_its+1;
    
    T_lambda = f(0) + epsilon*lambda*f_prime_0;
    T_alpha_lambda = f(0) + epsilon*alpha*lambda*f_prime_0;
  
    if isnan(f(lambda))
        lambda = 0
        s = 'f(lambda) is not a number';
        break
    end
      
    %vpa((f(h+lambda)-f(lambda))/(h),prec)
    
    if abs(vpa((f(h+lambda)-f(lambda))/(h),prec)) < 1e-22
        s = 'really small derivative';
        break
    end
    
    if lambda < 1e-22
        s = 'Lambda is really small';
        lambda = 0;
        break
    end
  
    num = 0;
    
    if f(lambda) > T_lambda
        lambda = lambda/alpha;
        num = num+1;
    end
    
    if f(lambda*alpha) < T_alpha_lambda 
        lambda = lambda*alpha;
        num = num+1;
    end
    
    if num == 0
        %s = 'ENDING'
        break
    end
end

if isnan(f(lambda)) || f(lambda)>f(0)
error('Bad job of the line search!')
end

end