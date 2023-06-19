function [lt,x] = Milstein_method(a,b,diff_b,T,N,x0,m1,m2)
h=T/N;
lt=0:h:T;
lx=zeros(size(lt));
lx(1)=x0;
Z=normal_generator(N,m1,m2);

for j=1:N
  lx(j+1) = lx(j) + h*a(lt(j),lx(j)) + b(lt(j),lx(j))*sqrt(h)*Z(j) + 0.5*b(lt(j),lx(j)) *diff_b(lt(j),lx(j))* ((sqrt(h)*Z(j))^2-h);
end
x = lx;
end