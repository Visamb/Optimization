%% nonlinearmin

x0 = [2;2];

[min,n] = nonlinearmin(@func,x0,1e-8,'BFGS')


%%

function [min, nbr_of_iterations] = nonlinearmin(f,x0,tol,method)

x_old = x0;
y = x_old;
dim = length(x_old);
D = eye(dim);
criterion = false;
y_new = 0
nbr_of_iterations = 0;

while criterion == false
nbr_of_iterations = nbr_of_iterations+1;

for j = 1:10
delta_f = grad(@func,y);
dj = -D*delta_f
F = @(lambda) func(y + dj*lambda);
[lambda,nbr,history] = linesearch(F,tol);
y_new = round(y + lambda*dj);
pj = lambda*dj;
qj = grad(@func,y_new)-grad(@func,y);

if method == 'BFGS'
D = D + (1+inv(transpose(pj)*qj)*transpose(qj)*D*qj)*inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(pj)*qj)*(pj*transpose(qj)*D + D*qj*transpose(pj));
else
D = D + inv(transpose(pj)*qj)*pj*transpose(pj)-inv(transpose(qj)*D*qj)*D*qj*transpose(qj)*D;
end
y = y_new
end

x_old = y;

if norm(grad(@func,x_old)) < tol
    criterion = true
end

end

min = x_old

end


