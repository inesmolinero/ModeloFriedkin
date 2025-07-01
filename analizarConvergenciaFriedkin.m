% analizarConvergenciaFriedkin.m
function out = analizarConvergenciaFriedkin(W,lambdas)
    n=numel(lambdas); A=(eye(n)-diag(lambdas))+diag(lambdas)*W;
    out.rhoA=max(abs(eig(A)));
    out.asimetria=norm(W-W','fro')/norm(W,'fro');
    v=eigs(W',1,'largestreal','Tolerance',1e-10); out.pi=abs(v)/sum(abs(v));
    M=W; prim=false; for k=1:50, M=M*W; if all(M(:)>0), prim=true; break; end; end; out.primitiva=prim;
end