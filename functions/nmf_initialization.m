function [W, H] = nmf_initialization(V, beta)
%nmf_initialization Summary of this function goes here
%   Detailed explanation goes here

s1 = size(V, 1);
s2 = size(V, 2);

W = rand(s1, 1);
H = rand(1, s2);

Wdiff = 0;
Hdiff = 0;

while Hdiff > 10^-7 && Wdiff > 10^7
    H1 = H;
    W1 = W;

    H = H.*(W'*((WH).^(beta-2).*V))./(W'*(WH).^(beta-1));
    W = W.*(((W*H).^(beta-2).*V)*H')./((W*H).^(beta-1)*H');
    
    Hdiff = norm(H1-H);
    Wdiff = norm(W1-W);
end
end

