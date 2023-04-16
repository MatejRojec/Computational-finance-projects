r = 0.06;
sigma = 0.3;
K = 10; 
a = 0;
S_star = 15;
T = 1;
NS = 300;
Nt = 10000;
hS = (S_star-a)/NS;
ht = T/Nt;


lS = linspace(a,S_star,NS);
lt = linspace(0,T,Nt);
grid = meshgrid(lS,lt);

u0 = @(S)(max(S-K,0));
ua = 0;
ub = @(t)(S_star - K*exp(-r*t));


grid(1,:) = u0(lS);
grid(:,1) = ua;
grid(:,end) = ub(lt);


i = 2:NS-1;
for j=1:Nt-1
        grid(j+1, i) = ht/2 * ((sigma^2) * i.^2 + r*i) .* grid(j,i+1) + (1-sigma^2 * i.^2 * ht - r*ht) .* grid(j,i) + ht/2 * (sigma^2 * i.^2 - r*i) .* grid(j,i-1);
end

% Create grid
[S, t] = meshgrid(lS, lt);
t = 1 - t;

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
title('European call option, V(S,t) explicit solutions')


figure(2)
surf(S,t,call_price, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option exact, V(S,t)')

figure(3)
surf(S,t,error, 'LineStyle','none', 'FaceColor','flat')
xlabel('S')
ylabel('t')
title('European call option, V(S,t) explicit method error')
max(max(error))