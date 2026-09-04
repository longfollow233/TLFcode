function pY = RBM(Xl,Yl,Xu,k)
U = size(Xu,1);
C = 2.^(-8:7);
lambda = 2.^(-8:7);
if U <= 1000
    mu = 2.^(-10:5);
else
    mu = ideal_sigma_2017(Xl, Yl);
end
num_C = length(C);
num_mu = length(mu);
num_lambda = length(lambda);
pY = zeros(U,num_C,num_mu,num_lambda);
ind = crossvalind('Kfold',U,k);
for p1 = 1:num_C
    for p2 = 1:num_mu
        H = zeros(U,U); 
        for i = 1:k
            teID = find(ind==i);
            trID = find(ind~=i);
            XXrr = kerf(Xu(trID,:),Xu(trID,:),mu(p2));
            XXre = kerf(Xu(trID,:),Xu(teID,:),mu(p2));
            Z = [XXrr+eye(length(trID))/C(p1),ones(length(trID),1); ones(1,length(trID)),0]\[XXre;ones(1,length(teID))];
            for i1 = 1:length(teID)
                for i2 = 1:length(trID)
                    H(teID(i1),trID(i2)) = H(teID(i1),trID(i2))+Z(i2,i1)/2;
                    H(trID(i2),teID(i1)) = H(trID(i2),teID(i1))+Z(i2,i1)/2;
                end
            end
        end
        mt = size(Xl,1);
        K = kerf(Xu,Xu,mu(p2));
        Ktt = kerf(Xu,Xl,mu(p2));
        g_base = [K+eye(U)/C(p1),ones(U,1);ones(1,U),0]\([Ktt;ones(1,mt)]*Yl);
        g_base = g_base(1:U, :);
        for p3 = 1:num_lambda
            g = g_base*lambda(p3);
            pY(:,p1,p2,p3) = IQPgreedyProMax(H,g); %>>>>>
        end
    end
end
end