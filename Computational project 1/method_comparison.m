r = 0.08;
sigma = 0.4;
K = 5; 
S_star = 40;
S = 5;
T = 4;

NS = [20 40 80 160 320];
Nt = [4000 8000 16000 32000 64000];

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
error_explicit_1 = abs(Values_1_explicit - exact_value_1)
Values_1_CN
error_CN_1 = abs(Values_1_CN - exact_value_1)





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
Error_2_explicit = abs(Values_2_explicit - exact_value_2)
Values_2_CN
Error_2_CN = abs(Values_2_CN - exact_value_2)


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
Error_3_explicit = abs(Values_3_explicit - exact_value_3)
Values_3_CN
Error_3_CN = abs(Values_3_CN - exact_value_3)


