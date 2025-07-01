% generarModeloFriedkin.m
% [W,x0,trolls,normales] = generarModeloFriedkin(p,sim,n,nneg,npos,or,seed)
function [W,x0,trolls,normales] = generarModeloFriedkin(p,sim,n,nneg,npos,op_range,seed)
    rng(seed,'twister');
    trolls = sort(randperm(n,nneg+npos));
    normales = setdiff(1:n,trolls);
    W = generarRedControlada(n,sim,p);
    x0 = zeros(n,1);
    x0(normales) = op_range(1)+diff(op_range)*rand(numel(normales),1);
    x0(trolls(1:nneg))=-1; x0(trolls(nneg+1:end))=1;
end