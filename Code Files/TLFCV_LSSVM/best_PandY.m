function res = best_PandY(Xl,Yl,Xu,Yu)
k = 5; % k = [5,10,15,20,LOO],select.
pY = RBM(Xl,Yl,Xu,k); %>>>>> 
[U,num_C,num_mu,num_lambda] = size(pY);
AC = zeros(num_C,num_mu,num_lambda);
for p1 = 1:num_C
    for p2 = 1:num_mu
        for p3 = 1:num_lambda
            AC(p1,p2,p3) = sum(pY(:,p1,p2,p3)==Yu) / U;
        end
    end
end
[max_value,max_index] = max(AC(:));
[best_C,best_mu,best_lambda] = ind2sub(size(AC),max_index);
res.py = pY(:,best_C,best_mu,best_lambda);
res.ac = max_value;
res.best_C = best_C-9;
res.best_mu = best_mu-11;
res.best_lambda = best_lambda-9;
res.AC = AC;
res.pY = pY;
end