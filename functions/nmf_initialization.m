function [W, H] = nmf_initialization(V, betaparam, stop, order)
%NMF_INITIALIZATION Summary of this function goes here
%   V : input matrix 
%   betaparam : as it is called in the algorithm
%   stop : stop parameter
%   OUTPUT W and H matrix

s1 = size(V, 1);
s2 = size(V, 2);

W = rand(s1, order);
H = rand(order, s2);

Vdiff = 1;
compteur=1;
while (Vdiff > stop) || compteur>1000
    compteur=compteur+1;
    
    H = H.*(W'*((W*H).^(betaparam-2).*V))./(W'*(W*H).^(betaparam-1));
    W = W.*(((W*H).^(betaparam-2).*V)*H')./((W*H).^(betaparam-1)*H');
    
    V2 = W*H;
    Vdiff = norm(V-V2);
end
end

