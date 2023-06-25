clear;
clc;

T=1;
N=1000;
x0 = 10;
h=T/N;
mu = 0.5;
sigma = 0.3;
a = @(t,x) mu * x;
b = @(t,x) sigma * x;
diff_b = @(t,x) sigma;
m1 = 0.2;
m2 = 0.5;

[lt,X] = Euler_Maruyama_method(a,b,T,N,x0,m1,m2);
[lt,X2] = Milstein_method(a,b,diff_b,T,N,x0,m1,m2);

Z=normal_generator(N,m1,m2);
S(1)=x0;
B = [0; cumsum(Z')] *sqrt(h);
for j=1:N
  S(j+1) = x0 * exp((mu - (sigma^2)/2) * lt(j+1) + sigma*B(j+1));
end

figure(1)
plot(lt,X,'b')
hold on
plot(lt,X2,'k')
plot(lt,S,'r')
legend({'Euler-Maruyama approximation','Milstein approximation','geometric Brownian motion'},'Location','best')
title('Comparison of theoretical solution and approximations')
xlabel('t')
ylabel('X(t)')
