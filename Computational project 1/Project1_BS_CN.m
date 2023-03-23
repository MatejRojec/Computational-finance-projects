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


u0 = @(S)(max(S-10,0));
ua = 0;
ub = @(t)(15 - 10*exp(-r*t));


grid(1,:) = u0(lS);
grid(:,1) = ua;
grid(:,end) = ub(lt);

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
B(2:NS-1:end-1) = -c(1:end-1);
B(NS-1:NS-1:end) = -a(2:end);
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
    u_solution = A \(B * vector + add);
    grid(j+1, 2:end-1) = u_solution;
end

lT = 1+zeros(1,Nt);
figure(1)
surf(lS,lT-lt,grid, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option, V(S,t)')
