function [g, l] = minimize(Z, t, w)
    
    % Av = b ...
    % v : concat(g, l);
    % size(g) = 256; size(l) = hw;

    lambda = 1;
    n = 256; % size of g vector
    [hw, k] = size(Z);
    
    A = sparse(hw*k + n + 1, n + hw);
    b = zeros(hw*k + n + 1, 1);

    % first term
    idx = 1;
    for i = 1:hw
        for j = 1:k
            weight = w(Z(i,j)+1);
            A(idx, Z(i,j)+1) = weight; 
            A(idx, n+i) = -weight;
            b(idx, 1) = weight * t(j);
            idx = idx + 1;
        end
    end

    % second term
    for i = 2:n-1
        A(idx, i-1) = lambda * w(i+1);
        A(idx, i) = -2 * lambda * w(i+1);
        A(idx, i+1) = lambda * w(i+1);
        idx = idx + 1;
    end

    v = A\b;
    g = v(1:n);
    l = v(n+1:end);
end