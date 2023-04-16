function [value] = exact(r, sigma, T, K, type, t, S)
% The function calculates the value of an option with the given paramaters
% obtained by the exact solutions 
% Inputs:
%   r - interest rate
%   sigma - volitility
%   T - expiration time
%   K - strike price
%   type - type of option either put or call
%   t - the time at witch we want to calculate the option price
%   S - the stock price at witch we want to calculate the stock price
% Outputs:
%   value - values of the option 

d1 = (log(S/K) + (r + 0.5*sigma^2)*(T-t)) ./ (sigma.*sqrt(T-t));
d2 = d1 - sigma.*sqrt(T-t);
Nd1 = normcdf(d1);
Nd2 = normcdf(d2);
Nmd1 = normcdf(-d1);
Nmd2 = normcdf(-d2);

if type == "call"
    value  = S.*Nd1 - K*exp(-r*(T-t)).*Nd2;
end

if type == "put"
    value = K*exp(-r*(T-t)).*Nmd2 - S.*Nmd1;
end

value = round(value,2);

end
