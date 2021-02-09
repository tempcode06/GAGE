function [U] = GAGE_EVD(Y1,Y2,R,thresh,iter)
%Geometry preserving Attributed Graph Embeddings
I = size(Y1,1);
Q0 = randn(I,R);
[Q0,~,~]=svd(Q0,'econ');
vec_ones=ones(I,1);
for i = 1:iter
    tmp1 = Q0 - 1/I*vec_ones*(vec_ones'*Q0);
    tmp1 = Y1*(Y1'*tmp1);
    tmp1 = tmp1 - 1/I*vec_ones*(vec_ones'*tmp1);
    tmp1 = Y1*(Y1'*tmp1);
    tmp1 = tmp1 - 1/I*vec_ones*(vec_ones'*tmp1);
    
    tmp2 = Q0 - 1/I*vec_ones*(vec_ones'*Q0);
    tmp2 = Y2*(Y2'*tmp2);
    tmp2 = tmp2 - 1/I*vec_ones*(vec_ones'*tmp2);
    tmp2 = Y2*(Y2'*tmp2);
    tmp2 = tmp2 - 1/I*vec_ones*(vec_ones'*tmp2);
    
    Z = tmp1 + tmp2;
    [Q,~] = qr(Z,0);
    if norm(Q-Q0,'fro')/norm(Q0,'fro')< thresh
        fprintf('convergence achieved in %i iterations \n', i)
        break
    end
    Q0 = Q;
end
tmp0 = Q - 1/I*vec_ones*(vec_ones'*Q);
tmp1 = Y1*(Y1'*tmp0);
tmp1 = tmp1 - 1/I*vec_ones*(vec_ones'*tmp1);
YY1 = Q'*tmp1;

tmp1 = Y2*(Y2'*tmp0);
tmp1 = tmp1 - 1/I*vec_ones*(vec_ones'*tmp1);
YY2 = Q'*tmp1;
% Retrieve U{1} using a GEVD.
[Apt,~] = eig(YY1.',YY2.');
U = cell(1,3);
U{1} = Q/(Apt.');
U{2}=U{1};

tmp = zeros(R,2);
for r=1:R
    tmp0 = U{1}(:,r) - 1/I*vec_ones*(vec_ones'*U{1}(:,r));
    tmp(r,1) = tmp0'*(Y1*(Y1'*tmp0));
    tmp(r,2) = tmp0'*(Y2*(Y2'*tmp0));
end

    
U{3} = (pinv((U{1}'*U{1}).^2)*tmp)';

vec_norm = cellfun(@(u) sqrt(sum(abs(u).^2)), U, 'UniformOutput', false);
U = cellfun(@(u,n) bsxfun(@rdivide, u, n), U, vec_norm, 'UniformOutput',0);
vec_norm = (prod(cat(1, vec_norm{:}),1)).^(1/3);
U = cellfun(@(u) bsxfun(@times, u, vec_norm), U, 'UniformOutput',0); 
end

