function [U] = GAGE(Y1,Y2,U,thresh,iter)
%Geometry preserving Attributed Graph Embeddings
I = size(Y1,1);
R = size(U{2},2);
vec_ones=ones(I,1);
for i=1:iter
    Uprev = U;
    tmp0 = (bsxfun(@times,U{2},U{3}(1,:)))';
    tmp1 = (bsxfun(@times,U{2},U{3}(2,:)))';
    tmp0 = tmp0 - 1/I*(tmp0*vec_ones)*vec_ones';
    tmp1 = tmp1 - 1/I*(tmp1*vec_ones)*vec_ones';
    tmp0 = ((tmp0*Y1)*Y1');
    tmp1 = ((tmp1*Y2)*Y2');
    tmp0 = tmp0 - 1/I*(tmp0*vec_ones)*vec_ones';
    tmp1 = tmp1 - 1/I*(tmp1*vec_ones)*vec_ones';
    tmp = tmp0 + tmp1;
    
    U{1} = (pinv((U{3}'*U{3}).*(U{2}'*U{2}))*tmp)';
    clear tmp tmp0 tmp1
    
    tmp0 = (bsxfun(@times,U{1},U{3}(1,:)))';
    tmp1 = (bsxfun(@times,U{1},U{3}(2,:)))';
    tmp0 = tmp0 - 1/I*(tmp0*vec_ones)*vec_ones';
    tmp1 = tmp1 - 1/I*(tmp1*vec_ones)*vec_ones';
    tmp0 = ((tmp0*Y1)*Y1');
    tmp1 = ((tmp1*Y2)*Y2');
    tmp0 = tmp0 - 1/I*(tmp0*vec_ones)*vec_ones';
    tmp1 = tmp1 - 1/I*(tmp1*vec_ones)*vec_ones';
    tmp = tmp0 + tmp1;
    
    U{2} = (pinv((U{3}'*U{3}).*(U{1}'*U{1}))*tmp)';
    clear tmp tmp0 tmp1

    tmp = zeros(R,2);
    for r=1:R
        tmp0 = U{1}(:,r) - 1/I*vec_ones*(vec_ones'*U{1}(:,r));
        tmp1 = U{2}(:,r) - 1/I*vec_ones*(vec_ones'*U{2}(:,r));
        tmp(r,1) = tmp0'*(Y1*(Y1'*tmp1));
        tmp(r,2) = tmp0'*(Y2*(Y2'*tmp1));
    end
    
    U{3} = (pinv((U{2}'*U{2}).*(U{1}'*U{1}))*tmp)';
    clear tmp tmp0 tmp1
    
    crit = [norm(Uprev{1}-U{1},'fro')/norm(Uprev{1},'fro'),norm(Uprev{2}-U{2},'fro')/norm(Uprev{2},'fro'),norm(Uprev{3}-U{3},'fro')/norm(Uprev{3},'fro')];
    
    if max(crit)<thresh
        fprintf('convergence achieved in %i iterations \n', i)
        break
    end
    
end
end

