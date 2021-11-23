%% nonlinearmin

x0 = [2;1];

%%

[x,no_its,normgrad] = nonlinearmin(@func,x0,1e-6,1);


%%

function [x, no_its, normg] = nonlinearmin(f,x0,tol,method)

%Stop criterion
criterion = false;

%Itilialize D
dim = length(x0);
D = eye(dim);

%Choose number of iterations for inner loop
n = length(x0);

%Initialize number of iterations for outer loop 
no_its = 0;

%First y (for inner loop)
y = x0;

%While stop criterion is not fulfilled
while criterion == false
no_its = no_its + 1


    for j = 1:n
        
        
        %Gradient at step y
        delta_f = grad(f,y)
        
        D = eye(dim)
        
        %Direction d
        dj = -D*delta_f
        
        %Add this to not get exploding hessians or impossible line search
        if norm(dj) < 10e-8
            break
        end
        
        %Linesearch to find optimal lambda
        F = @(lambda) f(y + dj*lambda);
        [lambda,nbr,history] = linesearch(F,1e-20);
        
        %Find new proposal for y
        y_new = (y + lambda*dj);
        
        %Construct pj and qj 
        pj = lambda*dj;
        qj = grad(f,y_new)-grad(f,y);

        %Update D with BFGS or DFP depending on chosen method (page 82,89)
        if method == 1
             D = D + (1+inv(transpose(pj)*qj)*transpose(qj)*D*qj)*inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(pj)*qj)*(pj*transpose(qj)*D + D*qj*transpose(pj))
        elseif method == 0
             D = D + inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(qj)*D*qj)*D*qj*transpose(qj)*D;
        end
        
        %Update y
        y = y_new;
        
        %Print 
        %fprintf('%12.4f %12.4f ',  norm(grad(@func,y)), lambda)


        
    end
    %When after inner loop, update x
    x = y;

if norm(grad(f,x)) < tol
    criterion = true;
end

end

normg = grad(f,x);

end


