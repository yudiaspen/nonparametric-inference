function jac = numeric_jacob(f, x0)
% Calculate Jacobian of function f at given x
epsilon = 1e-8; 
nx = length(x0); % Dimension of the input x;
f0 = f(x0); % caclulate f0, when no perturbation happens
jac=zeros(length(f0),nx);
for i = 1 : nx
    x = x0;
    x(i) =  x0(i)*(1+epsilon);
    jac(:, i) = (f(x) - f0) .* (1/epsilon/x0(i));
end
end