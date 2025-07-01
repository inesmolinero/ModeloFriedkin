function Wmod = ajustarTopologia(W,tipo)
    n=size(W,1);
    switch tipo
        case 'primitiva',    A=ones(n)+0.5*eye(n);
        case 'no_primitiva', A=diag(ones(n-1,1),1); A(n,1)=1;
        case 'conexa',       A=triu(ones(n),1)+triu(ones(n),1)';
        case 'no_conexa',    m=floor(n/2); A=blkdiag(ones(m)-eye(m),ones(n-m)-eye(n-m));
        otherwise, error('topología desconocida');
    end
    Wmod = normalizeRows(A);
end
