function res = best_PandY(Xl,Yl,Xu,Yu)
pY = BM(Xl,Yl,Xu); %%%%%>>>>>%%%%%>>>>>
[U,num_C,num_mu] = size(pY);
AC = zeros(num_C,num_mu);
for p1 = 1:num_C
    for p2 = 1:num_mu
        AC(p1,p2) = sum(pY(:,p1,p2) == Yu)/U;
    end
end
[best_C,best_mu] = find(AC == max(AC(:)));
if length(best_C)>1
    for h = 1:length(best_C)
        py = zeros(U,length(best_C));
        py(:,h) = pY(:,best_C(h),best_mu(h));
        res.py = py;
    end
else
    res.py = pY(:,best_C,best_mu);
end
res.ac = max(AC(:));
res.best_mu = best_mu-11;
res.best_C = best_C-9;
res.AC = AC;
res.pY = pY;
end