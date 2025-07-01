% generarRedControlada.m
function W = generarRedControlada(n,sim,p)
    A = rand(n)<p; A(1:n+1:end)=0;
    if sim, A = triu(A)+triu(A,1)'; end
    for i=1:n, if ~any(A(:,i)), j=randi(n); A(j,i)=1; if sim, A(i,j)=1; end; end; end
    W = normalizeRows(A);
end