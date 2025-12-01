

function [is_cons, is_pol, is_frag] = classify_state(x, labels, K, ...
        tau_cons, tau_gap, tau_size)
% Reglas exclusivas: consenso -> polarización (binaria genuina) -> fragmentación.
    if nargin<4, tau_cons = 1e-3; end
    if nargin<5, tau_gap  = 0.25; end
    if nargin<6, tau_size = 0.10; end

    % 1) Consenso (por rango)
    is_cons = (range(x) <= tau_cons) || (K==1);
    if is_cons
        is_pol = false; is_frag = false; return;
    end

    % 2) Polarización (exactamente 2 clústeres, bien separados y con tamaños comparables)
    if K == 2
        mu = accumarray(labels(:), x(:), [], @mean);
        nk = accumarray(labels(:), 1);
        gap = abs(mu(1)-mu(2));
        balance = min(nk)/numel(x);
        opposed = (sign(mu(1))*sign(mu(2)) == -1);
        is_pol = (gap >= tau_gap) && (balance >= tau_size) && opposed;
        is_frag = ~is_pol;      % exclusividad
        return;
    end

    % 3) Fragmentación (3 o más clústeres, o 2 que no cumplen polarización)
    is_pol = false;
    is_frag = (K >= 3);
end
