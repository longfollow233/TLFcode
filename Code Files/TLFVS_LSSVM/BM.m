function pY = BM(Xl,Yl,Xu)
L = size(Xl,1);
U = size(Xu,1);
C = 2.^(-8:7);
mu = 2.^(-10:5);
% mu = ideal_sigma_2017(Xl, Yl);
num_C  = length(C);
num_mu = length(mu);
pY = zeros(U,num_C,num_mu);
for p1 = 1:num_C
    for p2 = 1:num_mu
        K = kerf(Xu,Xu,mu(p2));
        Ktt = kerf(Xu,Xl,mu(p2));
        fval = [K+1/C(p1)*eye(U),ones(U,1);ones(1,U),0]\([Ktt;ones(1,L)]*Yl);
        pY(:,p1,p2) = sign(fval(1:U));
    end
end
end