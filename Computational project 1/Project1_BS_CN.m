% Crank-Nicolson method for European call option
clear all;
r = 0.06;
sigma = 0.3;
K = 10; 
a = 0;
S_star = 15;
T = 1;
NS = 50;
Nt = 10000;
hS = (S_star-a)/NS;
ht = T/Nt;

lS = linspace(a,S_star,NS);
lt = linspace(0,T,Nt);
grid = meshgrid(lS,lt);

% Initial conditions
u0 = @(S)(max(S-K,0));
ua = 0;
ub = @(t)(15 - 10*exp(-r*t));
grid(1,:) = u0(lS);
grid(:,1) = ua;
grid(:,end) = ub(lt);

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

% Create grid
[S, t] = meshgrid(lS, lt);
t = T - t;

% Compute Black-Scholes formula
d1 = (log(S/K) + (r + 0.5*sigma^2)*(T-t)) ./ (sigma.*sqrt(T-t));
d2 = d1 - sigma.*sqrt(T-t);
Nd1 = normcdf(d1);
Nd2 = normcdf(d2);
Nmd1 = normcdf(-d1);
Nmd2 = normcdf(-d2);
call_price = S.*Nd1 - K*exp(-r*(T-t)).*Nd2;
put_price = K*exp(-r*(T-t)).*Nmd2 - S.*Nmd1;

error = abs(call_price - grid);
lT = 1+zeros(1,Nt);
figure(1)
surf(lS,lT-lt,grid, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option, V(S,t)')

figure(2)
surf(S,t,call_price, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option exact, V(S,t)')

figure(3)
surf(S,t,error, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option CN, V(S,t) error')
max(max(error))