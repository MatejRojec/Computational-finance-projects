function [u] = uniform_generator(N, m1)
    M = 2^(31)-1;
    a = 16807;
    b = 0;
    m = zeros(1,N);
    m(1) = m1;
    for i=2:N
        m(i) = mod(a*m(i-1)+b,M);
    end
    u = m/M;
end
