function [vector] = ispec(X, window_length)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

vector = zeros(size(X,2)*window_length/2,1);

w = zeros(window_length, 1);
for k=1:window_length
    w(k) = sin((pi * k)/(window_length - 1)); % Sine-window calculation
end

for k = 1:size(X, 2)-1
    
    A = ifft(X(:, k));
    A = A.*w;
    vector((k-1)*window_length/2+1:(k-1)*window_length/2+window_length)=...
        vector((k-1)*window_length/2+1:(k-1)*window_length/2+window_length)+A;
    
end



end