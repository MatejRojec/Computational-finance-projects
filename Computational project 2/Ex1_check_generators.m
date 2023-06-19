N=1000000;
U = uniform_generator(N,0.2);

figure(1)
[g1,x1] = ksdensity(U);
f1 = unifpdf(x1,0,1);
plot(x1,f1,'-r','LineWidth',4)
hold on
histogram(U,'Normalization','pdf')
legend({'theoretical pdf','histogram'},'Location','best')
title('Theoretical probability density function and empirical histogram')
xlabel('x')
ylabel('f(x)')



X = normal_generator(N,0.2,0.5);
figure(2)
[g2,x2] = ksdensity(X);
f2 = normpdf(x2,0,1);
plot(x2,f2,'-r','LineWidth',4)
hold on
histogram(X,'Normalization','pdf')
legend({'theoretical pdf','histogram'},'Location','best')
title('Theoretical probability density function and empirical histogram')
xlabel('x')
ylabel('f(x)')