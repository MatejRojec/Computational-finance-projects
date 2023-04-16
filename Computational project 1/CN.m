function [value] = CN(r, sigma, T, K, S_star, NS, Nt, type, time, stock)
% The function calculates the value of an option with the given paramaters
% obtained by the CN method 
% Inputs:
%   r - interest rate
%   sigma - volitility
%   T - expiration time
%   K - strike price
%   S_star - boundary condition limit
%   NS - number of points on the S axis
%   Nt - number of points on the t axis
%   type - type of option either put or call
%   time - the time at witch we want to calculate the option price
%   stock - the stock price at witch we want to calculate the stock price
% Outputs:
%   value - values of the option 

a = 0;
hS = (S_star-a)/NS;
ht = T/Nt;
lS = linspace(a,S_star,NS);
lt = linspace(0,T,Nt);
grid = meshgrid(lS,lt);

% Initial conditions
if type == "call"
    u0 = @(S)(max(S-K,0));
    ua = 0;
    ub = @(t)(S_star - K*exp(-r*t));
    grid(1,:) = u0(lS);
    grid(:,1) = ua;
    grid(:,end) = ub(lt);
end

if type == "put"
    u0 = @(S)(max(K-S,0));
    ua = 0;
    ub = @(t)(K*exp(-r*t));
    grid(1,:) = u0(lS);
    grid(:,1) = ua;
    grid(:,end) = ub(lt);
end

% Crank-Nicolson method
lambda = ht/hS^2;
i = 1:NS-2;
a = -(sigma^2)/4 * ht * i.^2 + r*ht/4*i;
b = 1 +(sigma^2)/2 * ht * i.^2 + r*ht/2;
c = -(sigma^2)/4 * ht * i.^2 - r*ht/4*i;
d = 1 -(sigma^2)/2 * ht * i.^2 - r*ht/2;
A = zeros(NS-2);
A(1:NS-1:end) = b;
A(2:NS-1:end-1) = a(2:end);
A(NS-1:NS-1:end) = c(1:end-1);
B = zeros(NS-2);
B(1:NS-1:end) = d;
B(2:NS-1:end-1) = -a(2:end);
B(NS-1:NS-1:end) = -c(1:end-1);
first_col = zeros(1,NS-2);
last_col = zeros(1,NS-2);
first_col(1) = -a(1);
last_col(end) = -c(end);
B = [first_col' B last_col'];
for j=1:Nt-1
    vector = grid(j, :)';
    add = zeros(1,NS-2)';
    add(1) = -a(1)* grid(j+1,1);
    add(end) = -c(end)*grid(j+1, end);

    % add boundary conditions
    vector(1) = ua;
    vector(end) = ub(lt(j+1));
    u_solution = A \(B * vector + add);    
    % insert boundary conditions into solution
    u_solution = [ua; u_solution; ub(lt(j+1))];    
    grid(j+1, :) = u_solution;
end

idx_t = find(lt >= T - time, 1, 'first');
idx_s = find(lS >= stock, 1, 'first');

value = round(grid(idx_t, idx_s), 2);


end
