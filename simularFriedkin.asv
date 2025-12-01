function X = simularFriedkin(W, x0, lambdas, T, x_star, tol_fp)
    n = numel(x0);
    X = zeros(n, T);
    x = x0;
    Lam = spdiags(lambdas,0,n,n);
    has_target = (nargin >= 5) && ~isempty(x_star) && (nargin >= 6) && ~isempty(tol_fp);
    for t = 1:T
        x = Lam*W*x + (speye(n)-Lam)*x0;
        X(:,t) = x;
        if has_target
            if norm(x - x_star, Inf) < tol_fp
                X = X(:,1:t);   % recorta y termina
                return
            end
        end
    end
end
