% simularFriedkin.m
function x_hist = simularFriedkin(W,x0,lambdas,n_iters)
    n=numel(x0); x=x0; x_hist=zeros(n,n_iters);
    for t=1:n_iters
        x=(eye(n)-diag(lambdas))*x + diag(lambdas)*W*x;
        x_hist(:,t)=x;
    end
end