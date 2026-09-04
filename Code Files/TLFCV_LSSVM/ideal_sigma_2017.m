function mu = ideal_sigma_2017(X, Y)
% The best sigma of the Gaussian kernel is given according to the ideal kernel theory.

Mu = 2.^(-10:5);
num_X = size(X, 1);

J = double(Y == Y'); % 显式转换为数值矩阵 (增强健壮性)

num_mu = length(Mu);
Value = zeros(num_mu, 1);

D = pdist2(X, X, 'euclidean');
for p1 = 1:num_mu
    K = exp(-Mu(p1) * D.^2);
    Value(p1) = sum(sum(abs(K - J))) / (num_X^2); % 关键修改点：双重sum求和
end

[~, index] = min(Value);
mu = Mu(index);
end