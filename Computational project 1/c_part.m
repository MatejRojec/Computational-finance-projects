r = 0.08;
sigma = 0.4;
K = 5; 
S_star = 40;
S = 5;
T = 4;


NS = [10 20 40 80 160 300 1000];
Nt = [1000 2000 3000 4000 5000 10000];

time = 1;
exact_value_1 = exact(r, sigma, T, K, "put", time, S);
Values_1_explicit = zeros(length(NS), length(Nt));
Values_1_CN = zeros(length(NS), length(Nt));
for i =1:length(NS)
    for j = 1:length(Nt)
        Values_1_explicit(i, j) = explicit(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S);    
        Values_1_CN(i, j) = CN(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S); 
    end
end
Values_1_explicit
Values_1_CN
exact_value_1




time = 2;
exact_value_2 = exact(r, sigma, T, K, "put", time, S);
Values_2_explicit = zeros(length(NS), length(Nt));
Values_2_CN = zeros(length(NS), length(Nt));
for i =1:length(NS)
    for j = 1:length(Nt)
        Values_2_explicit(i, j) = explicit(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S);    
        Values_2_CN(i, j) = CN(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S); 
    end
end
Values_2_explicit
Values_2_CN
exact_value_2


time = 3;
exact_value_3 = exact(r, sigma, T, K, "put", time, S);
Values_3_explicit = zeros(length(NS), length(Nt));
Values_3_CN = zeros(length(NS), length(Nt));
for i =1:length(NS)
    for j = 1:length(Nt)
        Values_3_explicit(i, j) = explicit(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S);    
        Values_3_CN(i, j) = CN(r, sigma, T, K, S_star, NS(i), Nt(j), "put", time, S); 
    end
end
Values_3_explicit
Values_3_CN
exact_value_3


