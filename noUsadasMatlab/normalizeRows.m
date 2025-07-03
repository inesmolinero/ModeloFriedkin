% normalizeRows.m
function M = normalizeRows(A)
    s=sum(A,2); s(s==0)=1; M=A./s;
end
