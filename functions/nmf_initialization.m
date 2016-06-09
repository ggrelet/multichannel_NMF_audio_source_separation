function [W, H] = nmf_initialization(V, betaparam, stop)
%NMF_INITIALIZATION Summary of this function goes here
%   V : input matrix 
%   betaparam : as it is called in the algorithm
%   stop : stop parameter
%   OUTPUT W and H matrix

s1 = size(V, 1);
s2 = size(V, 2);

W = rand(s1, 1);
H = rand(1, s2);

Wdiff = 1;
Hdiff = 1;

while Hdiff > stop && Wdiff > stop
    H1 = H;
    W1 = W;

    H = H.*(W'*((W*H).^(betaparam-2).*V))./(W'*(W*H).^(betaparam-1));
    W = W.*(((W*H).^(betaparam-2).*V)*H')./((W*H).^(betaparam-1)*H');
    
    Hdiff = norm(H1-H);
    Wdiff = norm(W1-W);
end
end

