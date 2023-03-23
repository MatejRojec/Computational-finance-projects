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
[X,T] = meshgrid(lS,lt);

u0 = @(S)(max(S-10,0));
ua = 0;
ub = @(t)(15 - 10*exp(-r*t));


grid(1,:) = u0(lS);
grid(:,1) = ua;
grid(:,end) = ub(lt);


i = 2:NS-1;
for j=1:Nt-1
        grid(j+1, i) = ht/2 * ((sigma^2) * i.^2 + r*i) .* grid(j,i+1) + (1-sigma^2 * i.^2 * ht - r*ht) .* grid(j,i) + ht/2 * (sigma^2 * i.^2 - r*i) .* grid(j,i-1);
end


lT = 1+zeros(1,Nt);
figure(1)
surf(lS,lT-lt,grid, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option, V(S,t)')
