function x = IQPgreedyProMax(Q,f)
%IQPSOLVER consider integer quadratic programming:
%  max x'*Q*x+f'x
%  s.t. x= +1 or -1...
%  **where the diagonals of Q are 0, Q=Q'.

m=length(f);
x=zeros(m,1);
ms=m*(m+1)/2; % m+m(m-1)/2

% 1. build s and sort |s|
s=[zeros(1,ms-m),f'];
row=[zeros(1,ms-m),1:m];
col=zeros(1,ms);
t=1;
for i=1:m
    for j=i+1:m
        row(1,t)=i;
        col(1,t)=j;
        s(1,t)=Q(i,j)*2;
        t=t+1;
    end
end

% Remove tiny values to speed up
ind=find(s<1e-8 & s>-1e-8);
if ~isempty(ind)
    s(:,ind)=[];
    row(:,ind)=[];
    col(:,ind)=[];
    ms=ms-length(ind);
end

[~,ind]=sort(abs(s),'descend'); %ss=abs(s)
s=s(1,ind);
row=row(1,ind);
col=col(1,ind);

% 2. set solution
% graph, 1st is positive connect set, 2nd is negative connect set, 3rd is
%     bag-value, 4th is all of the pairs index, and 5th is the first contradiction point
G={};
bg=BagGraph(G);
% disp(class(bg))
t=1;

while(t<=ms)
    if col(1,t)~=0 % Case: Pairwise Term (h_ij)                    
        if x(row(1,t),1)==0 && x(col(1,t),1)==0 % both are no label
           isUpdateS=updateGraph(bg,row(1,t),col(1,t),s(1,t),t);
           if isUpdateS>0 % need goto tth step
               s(1,bg.G{isUpdateS,4})=0;
               t=bg.G{isUpdateS,5};
               bg.G(isUpdateS,:)=cell(1,5);
               continue;
           end
        elseif x(row(1,t),1)==0 % the other has a label
            if s(1,t)>0
                x(row(1,t),1)=x(col(1,t),1);
            elseif s(1,t)<0
                x(row(1,t),1)=-x(col(1,t),1);
            end
            % Propagate to bag
            i1=0;
            mG=size(bg.G,1);
            for i=1:mG
                if ~isempty(find(bg.G{i,1}==row(1,t),1))
                    i1=i; % in left bag
                    break;
                end
                if ~isempty(find(bg.G{i,2}==row(1,t),1))
                    i1=mG+i; % in right bag
                    break;
                end
            end
            if i1>mG
                if ~isempty(bg.G{i1-mG,1})
                    x(bg.G{i1-mG,1},1)=-x(row(1,t),1);
                end
                x(bg.G{i1-mG,2},1)=x(row(1,t),1);
                bg.G(i1-mG,:)=cell(1,5);
            elseif i1>0
                x(bg.G{i1,1},1)=x(row(1,t),1);
                if ~isempty(bg.G{i1,2})
                    x(bg.G{i1,2},1)=-x(row(1,t),1);
                end
                bg.G(i1,:)=cell(1,5);
            end
        elseif x(col(1,t),1)==0 % the other has a label
            if s(1,t)>0
                x(col(1,t),1)=x(row(1,t),1);
            elseif s(1,t)<0
                x(col(1,t),1)=-x(row(1,t),1);
            end
            % Propagate to bag
            i1=0;
            mG=size(bg.G,1);
            for i=1:mG
                if ~isempty(find(bg.G{i,1}==col(1,t),1))
                    i1=i; % in left bag
                    break;
                end
                if ~isempty(find(bg.G{i,2}==col(1,t),1))
                    i1=mG+i; % in right bag
                    break;
                end
            end
            if i1>mG
                if ~isempty(bg.G{i1-mG,1})
                    x(bg.G{i1-mG,1},1)=-x(col(1,t),1);
                end
                x(bg.G{i1-mG,2},1)=x(col(1,t),1);
                bg.G(i1-mG,:)=cell(1,5);
            elseif i1>0
                x(bg.G{i1,1},1)=x(col(1,t),1);
                if ~isempty(bg.G{i1,2})
                    x(bg.G{i1,2},1)=-x(col(1,t),1);
                end
                bg.G(i1,:)=cell(1,5);
            end
        end
        
    else    % Case: Single Term (g_i) - >>> 关键修改区域 <<<
        
        idx = row(1,t);
        % 【FIX START】: 只有当节点还没有被之前的强连接(H)赋值时，才使用 g 进行赋值
        if x(idx,1) == 0 
            if s(1,t)>0
                x(idx,1)=1;
            else
                x(idx,1)=-1;
            end
            
            % 只有发生了新赋值，才需要进行 Bag 传播
            i1=0;
            mG=size(bg.G,1);
            for i=1:mG
                if ~isempty(find(bg.G{i,1}==idx,1))
                    i1=i; % in left bag
                    break;
                end
                if ~isempty(find(bg.G{i,2}==idx,1))
                    i1=mG+i; % in right bag
                    break;
                end
            end
            if i1>mG
                if ~isempty(bg.G{i1-mG,1})
                    x(bg.G{i1-mG,1},1)=-x(idx,1);
                end
                x(bg.G{i1-mG,2},1)=x(idx,1);
                bg.G(i1-mG,:)=cell(1,5);
            elseif i1>0
                x(bg.G{i1,1},1)=x(idx,1);
                if ~isempty(bg.G{i1,2})
                    x(bg.G{i1,2},1)=-x(idx,1);
                end
                bg.G(i1,:)=cell(1,5);
            end
        end
        % 【FIX END】: 如果 x(idx) 已经有值了，说明它被更强的 H 项决定了，跳过弱 g 项
    end            
    t=t+1;
end

% step 3, fill the rest (clean up any remaining uncertain bags)
mG=size(bg.G,1);
for i=1:mG
    if ~isempty(bg.G{i,1})
        x(bg.G{i,1},1)=1;
    end
    if ~isempty(bg.G{i,2})
        x(bg.G{i,2},1)=-1;
    end
end
% Check if any variables remain 0 (should not happen usually, but good for safety)
unassigned = find(x==0);
if ~isempty(unassigned)
    % Assign remaining based on their individual g (original f)
    % Note: f is not sorted here, need original index if needed, 
    % but strictly speaking they should have been handled.
    % Just random or default to 1 if absolutely no info.
    x(unassigned) = 1; 
end

end

function isUpdateS=updateGraph(bg,y1,y2,st,t)
% *For two points have no label!
% bg is the bag-graph, y1 and y2 are a pair of indices, st is the weight of
% pair
isUpdateS=0;  % 0 -- no change, >0 -- the ith row that have to delete
if st==0
    return;
end
m=size(bg.G,1);
i1=0;
i2=0;
for i=1:m
    if ~isempty(find(bg.G{i,1}==y1,1))
        i1=i; % in left bag
        break;
    end
    if ~isempty(find(bg.G{i,2}==y1,1))
        i1=m+i; % in right bag
        break;
    end
end
for i=1:m
    if ~isempty(find(bg.G{i,1}==y2,1))
        i2=i; % in left bag
        break;
    end
    if ~isempty(find(bg.G{i,2}==y2,1))
        i2=m+i; % in right bag
        break;
    end
end
if i1==0 && i2==0 % (a) 1
    bg.G(end+1,:)=cell(1,5);
    bg.G{m+1,1}=[y1,y2];
    bg.G{m+1,3}=abs(st);
    bg.G{m+1,4}=t;
    bg.G{m+1,5}=0;
elseif i1==0 % (a) 2
    if st>0
        if i2>m
            bg.G{i2-m,2}(1,end+1)=y1;
            bg.G{i2-m,3}=bg.G{i2-m,3}+st;
            bg.G{i2-m,4}(1,end+1)=t;
        else
            bg.G{i2,1}(1,end+1)=y1;
            bg.G{i2,3}=bg.G{i2,3}+st;
            bg.G{i2,4}(1,end+1)=t;
        end 
    else 
        if i2>m
            bg.G{i2-m,1}(1,end+1)=y1;
            bg.G{i2-m,3}=bg.G{i2-m,3}-st;
            bg.G{i2-m,4}(1,end+1)=t;
        else
            bg.G{i2,1}(1,end+1)=y1;
            bg.G{i2,3}=bg.G{i2,3}-st;
            bg.G{i2,4}(1,end+1)=t;
        end
    end
elseif i2==0 % (a) 3
    if st>0
        if i1>m
            bg.G{i1-m,2}(1,end+1)=y2;
            bg.G{i1-m,3}=bg.G{i1-m,3}+st;
            bg.G{i1-m,4}(1,end+1)=t;
        else
            bg.G{i1,1}(1,end+1)=y2;
            bg.G{i1,3}=bg.G{i1,3}+st;
            bg.G{i1,4}(1,end+1)=t;
        end 
    else 
        if i1>m
            bg.G{i1-m,1}(1,end+1)=y2;
            bg.G{i1-m,3}=bg.G{i1-m,3}-st;
            bg.G{i1-m,4}(1,end+1)=t;
        else
            bg.G{i1,1}(1,end+1)=y2;
            bg.G{i1,3}=bg.G{i1,3}-st;
            bg.G{i1,4}(1,end+1)=t;
        end
    end
else
    if i1==i2 % in the same bag: 1 or 2
        if st>0 % no contradiction
            if i1>m
                bg.G{i1-m,3}=bg.G{i1-m,3}+st;
                bg.G{i1-m,4}(1,end+1)=t;
            else
                bg.G{i1,3}=bg.G{i1,3}+st;
                bg.G{i1,4}(1,end+1)=t;
            end
        else  % rises contradiction
            if i1>m
                bg.G{i1-m,3}=bg.G{i1-m,3}+st;
                if bg.G{i1-m,5}==0
                    bg.G{i1-m,5}=t;
                end
                if bg.G{i1-m,3}<0
                    isUpdateS=i1-m;
                end
            else
                bg.G{i1,3}=bg.G{i1,3}+st;
                if bg.G{i1,5}==0
                    bg.G{i1,5}=t;
                end
                if bg.G{i1,3}<0
                    isUpdateS=i1;
                end
            end
        end
    elseif i1-i2==m % in the same bag: i1--2, i2--1
        if st>0 % rise contradiction
            bg.G{i2,3}=bg.G{i2,3}-st;
            if bg.G{i2,5}==0
                bg.G{i2,5}=t;
            end
            if bg.G{i2,3}<0
                isUpdateS=i2;
            end
        else  % no contradiction
            bg.G{i2,3}=bg.G{i2,3}-st;
            bg.G{i2,4}(1,end+1)=t;
        end
    elseif i2-i1==m % in the same bag: i1--1, i2--2
        if st>0 % rise contradiction
            bg.G{i1,3}=bg.G{i1,3}-st;
            if bg.G{i1,5}==0
                bg.G{i1,5}=t;
            end
            if bg.G{i1,3}<0
                isUpdateS=i1;
            end
        else
            bg.G{i1,3}=bg.G{i1,3}-st;
            bg.G{i1,4}(1,end+1)=t;
        end
    else % in two different bags
        if i1>m && i2>m
            imin=min(i1,i2)-m;
            imax=max(i1,i2)-m;
        elseif i1>m
            imin=min(i1-m,i2);
            imax=max(i1-m,i2);
        elseif i2>m
            imin=min(i1,i2-m);
            imax=max(i1,i2-m);
        else
            imin=min(i1,i2);
            imax=max(i1,i2);
        end
        if (i1>m && i2>m) || (i1<=m && i2<=m)
            if st>0
                bg.G{imin,1}=[bg.G{imin,1},bg.G{imax,1}]; 
                bg.G{imin,2}=[bg.G{imin,2},bg.G{imax,2}]; 
                bg.G{imin,3}=bg.G{imin,3}+bg.G{imax,3}+st;
            else
                bg.G{imin,1}=[bg.G{imin,1},bg.G{imax,2}]; 
                bg.G{imin,2}=[bg.G{imin,2},bg.G{imax,1}]; 
                bg.G{imin,3}=bg.G{imin,3}+bg.G{imax,3}-st;
            end
        else
            if st>0
                bg.G{imin,1}=[bg.G{imin,1},bg.G{imax,2}]; 
                bg.G{imin,2}=[bg.G{imin,2},bg.G{imax,1}];
                bg.G{imin,3}=bg.G{imin,3}+bg.G{imax,3}+st;
            else
                bg.G{imin,1}=[bg.G{imin,1},bg.G{imax,1}]; 
                bg.G{imin,2}=[bg.G{imin,2},bg.G{imax,2}]; 
                bg.G{imin,3}=bg.G{imin,3}+bg.G{imax,3}-st;
            end                
        end
        bg.G{imin,4}=[bg.G{imin,4},bg.G{imax,4}]; 
        bg.G{imin,5}=min(bg.G{imin,5},bg.G{imax,5});
        bg.G(imax,:)=cell(1,5);
    end
end                                            
end