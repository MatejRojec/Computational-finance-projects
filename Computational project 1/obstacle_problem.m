% function g 
g = @(x) 15/16 + 3/8*x - 25/16*x.^2;
% distritise the interval 

N = 200; % number of grid points
x = linspace(-1,1,N)';
x=x(2:(end-1));
h = x(2) - x(1);

% Define the diagonal, upper and lower diagonals of the matrix M as follows
a = 2 * ones(N-2, 1);
b = -1 * ones(N-3,1);
A = diag(a) + diag(b, 1) + diag(b, -1)

% right hand side of b 
g_matrix = - A * g(x);

% initialise x0
x0 = zeros(N-2,1);

% parametrs 
omega = 1;
maxit = 8259;
n = 0

x = x0;
while n < maxit
    for i = 1:N-2
        mi = A(i,:);
        r = g_matrix(i) - mi*x;
        x(i) = max(0, x(i) + omega/A(i,i)* r);
    end
    n = n + 1;
end

% plot 

xx = linspace(-1,1,200)';
xx = xx(2:(end-1));

plot(xx', g(xx) ,'b', xx', x + g(xx), 'r');
xlim([-1,1]);
xlabel('x');
