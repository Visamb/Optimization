%% nonlinearmin

%%


%[x, no_its, normg] = nonlinearmin(@func,[-20;20],1e-6,1,0)

%%

function [x, no_its, normg] = nonlinearmin(f,x0,tol,method,restart)

%Stop criterion
criterion = false;

%Itilialize D
dim = length(x0);
D = eye(dim);

%Choose number of iterations for inner loop
n = length(x0);

%Initialize number of iterations for outer loop 
no_its = 0;
ls_its = 0;

%First y (for inner loop)
y = x0;



%Create table for printing
fileID = fopen('data.txt','w');
fprintf('%12s %12s %12s %12s %12s %12s %12s\n', 'iteration','x','step size', 'f(x)', 'norm(grad)', 'ls iters', 'lambda');
format short g
%While stop criterion is not fulfilled
while criterion == false

    
    no_its = no_its + 1;
   
    if restart
            D = eye(dim);
    end

    for j = 1:n
        
        
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
        [lambda,ls_its] = armijo(F,2,0.01);
        
        %Find new proposal for y
        y_new = (y + lambda*dj);
        
        %Construct pj and qj 
        pj = lambda*dj;
        qj = grad(f,y_new)-grad(f,y);

        %Update D with BFGS or DFP depending on chosen method (page 82,89)
        if method == 1
             D = D + (1+inv(transpose(pj)*qj)*transpose(qj)*D*qj)*inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(pj)*qj)*(pj*transpose(qj)*D + D*qj*transpose(pj));
        elseif method == 0
             D = D + inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(qj)*D*qj)*D*qj*transpose(qj)*D;
        end
        
        %Update y
        y_prev = y;
        y = y_new;
        
        grad_norm = norm(grad(f,y));
        step_size = abs(dj*lambda);
        
        %Print 
        %fprintf('%12.4f %12.4f ',  norm(grad(@func,y)), lambda)
        
    end
    
    %Print
    fprintf('%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', no_its,y(1)',norm(step_size),f(y),grad_norm,ls_its,lambda)
    
    for j = 1:dim-1
        fprintf('%12s %12.4f\n', strings,y(j+1))
    end


    %When after inner loop, update x
    x = y;
    grads = grad(f,x);

if norm(grad(f,x)) < tol || abs(f(y) - f(y_prev)) < 1e-6
    criterion = true;
end

end



normg = grad(f,x);



end


