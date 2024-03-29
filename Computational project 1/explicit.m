function [value] = explicit(r, sigma, T, K, S_star, NS, Nt, type, time, stock)
% The function calculates the value of an option with the given paramaters
% obtained by the explicit method 
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

i = 2:NS-1;
for j=1:Nt-1
        grid(j+1, i) = ht/2 * ((sigma^2) * i.^2 + r*i) .* grid(j,i+1) + (1-sigma^2 * i.^2 * ht - r*ht) .* grid(j,i) + ht/2 * (sigma^2 * i.^2 - r*i) .* grid(j,i-1);
end

idx_t = find(lt >= T - time, 1, 'first');
idx_s = find(lS >= stock, 1, 'first');

value = round(grid(idx_t, idx_s), 2);

end
