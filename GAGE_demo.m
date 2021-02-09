%Geometry preserving Attributed Graph Embeddings
clear all;clc;close all;
data = 'webkb'; % webkb, wiki, BlogCatalog
load(strcat('data/',strcat(data),'.mat'))

%% get embeddings from CP demcomposition
dim = 128; %embedding dimension
tic
iter1 = 200;
thresh1 = 0.01;
[Ugevd] = GAGE_EVD(network,sparse(feature),dim,thresh1,iter1);
iter2 = 200;
thresh2 = 0.01;
[U] = GAGE(network,sparse(feature),Ugevd,thresh2,iter2);
toc
lamda = ones(dim,1);
for i=1:3
    for j=1:dim
        tmp =norm(U{i}(:,j));
        U{i}(:,j) = U{i}(:,j)/tmp;
        lamda(j) = lamda(j)*tmp;
    end
end
C=U{3}*diag(lamda);
tmp2 = sign(C(1,:));
lambda = 0.8; % parameter that weights the contribution of connectivity
% attributed info
emb = U{1}*diag(tmp2)*diag(sqrt(lambda*abs(C(1,:))+(1-lambda)*abs(C(2,:)))); %embeddings


