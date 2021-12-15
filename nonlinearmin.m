%% nonlinearmin
x0 = [200,200]

[x0, no_its, normg] = nonlinearmin(@rosenbrock,x0,1e-6,1,0,1);


%%


function [x, no_its, grad_norm] = nonlinearmin(f,x0,tol,method,restart,printout)
% Optimization based on quasi-Newton optimization algorithms (BFGS or DFP).
%
% :param f: input function
% :param x0: vector, intial start point(s).
% :param tol: double, tolerance for termination criteria.
% :param restart: boolean, restart if true
% :param method: boolean, choice of algorithm, true for 'BFGS and false for 'DFP'.
% :param printout: integer/boolean, prints status if true.
% :return: x: vector of doubles, found optimal point(s).
% :return: no_its: integer, number of iterations.
% :return: grad_norm: double, gradient of the norm at termination.

%Stop criterion
criterion = false;

%Itilialize D
dim = length(x0);

%Choose number of iterations for inner loop
n = length(x0);

%Initialize number of iterations for outer loop 
no_its = 0;
ls_its = 0;

%Initialize parameters in case of immediate stop.
lambda = 0;
step_size = 0;
grad_norm = 0;

%First y (for inner loop)
y = x0;

%Create table for printing
if printout
    fprintf('%12s %12s %12s %12s %12s %12s %12s\n', 'iteration','x','step size', 'f(x)', 'norm(grad)', 'ls iters', 'lambda');
    fprintf('%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', no_its,y(1)',norm(step_size),f(y),norm(grad(f,y)),0,0);
    for l = 1:dim-1
        fprintf('%12s %12.4f\n', strings,y(l+1))
    end
    format short g
end

%While stop criterion is not fulfilled
    while criterion == false
        
        D = eye(dim);
        
        %In case of restart, start from gradient descent at some m<n.
        if restart
            n = ceil(n/1.8);    
        end
        
        ls_its = 0;

        %Inner loop
        for j = 1:n
            
            %Increase counter.
            no_its = no_its + 1;

            %Gradient at step y
            delta_f = grad(f,y);

            %Direction d
            dj = -D*delta_f;

            %Add this to not get exploding hessians or impossible line search
            if norm(dj) < 10e-8
                break
            end

            %Linesearch to find optimal lambda
            F = @(lamb) f(y + dj*lamb);
            f_prime_0 = delta_f'*dj;
            [lambda,ls_its] = armijo2(F,2,0.05,f_prime_0);

            %Find new proposal for y
            y_new = (y + lambda*dj);
            delta_f_new = grad(f,y_new);

            %Construct pj and qj 
            pj = lambda*dj;
            qj = delta_f_new-delta_f;

            %Update D with BFGS or DFP depending on chosen method (page 82,89)
            if method == 1
                 D = D + (1+(transpose(pj)*qj)^(-1)*transpose(qj)*D*qj)*(transpose(pj)*qj)^(-1)*pj*transpose(pj)-(transpose(pj)*qj)^(-1)*(pj*transpose(qj)*D + D*qj*transpose(pj));
            else 
                 D = D + (transpose(pj)*qj)^(-1)*pj*transpose(pj)-(transpose(qj)*D*qj)^(-1)*D*qj*transpose(qj)*D;
            end

            %Update y
            y_prev = y;
            y = y_new;

            %Norm of gradient and last step size.
            grad_norm = norm(delta_f_new);
            step_size = abs(pj);
            
             if printout
                %Print status
                fprintf('%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', no_its,y(1)',norm(step_size),f(y),grad_norm,ls_its,lambda)
                for l = 1:dim-1
                    fprintf('%12s %12.4f\n', strings,y(l+1))
                end
             end
             
                 %Check criterion
            if grad_norm < tol 
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                fprintf('%s\n', 'Optimization was terminated.')
                fprintf('%s, %f.\n ','Stationary point reached. Gradient was smaller than the specified tolerance', tol)
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                criterion = true;
                break
            end

             %Check criterion
            f_y_prev = f(y_prev);
            if  (abs(f(y) - f_y_prev)/abs(f(y_prev)) < tol*1e-9)
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                fprintf('%s\n', 'Optimization was terminated.')
                fprintf('%s, %f.\n ','Stopping criterion: Change in function value was less than the specified tolerance', tol)
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                criterion = true;
                break
            end

             %Check criterion
            if (norm(y-y_prev) < tol*1e-9)
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                fprintf('%s\n', 'Optimization was terminated.')
                fprintf('%s, %f.\n ','Stopping criterion: Step size was less than the specified tolerance', tol)
                fprintf('%s\n', '------------------------------------------------------------------------------------------')
                criterion = true;
                break
            end
        end

        %When after inner loop, update x
        x = y;

    end

end


