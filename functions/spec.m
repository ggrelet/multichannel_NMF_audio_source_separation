function [X] = spec(vector, window_length, overlap)
%SPEC Spectrogram function with sine-window
%   vector : input x vector
%   window_length : size of the sine-window, typically 1024 samples
%   overlap : overlapping coefficient, typically 50%
%   OUTPUT : 

vector_length = size(vector(:, 1), 1); % Size of the input vector (the sound)

w = zeros(window_length, 1); % Window vector
X = zeros(window_length, ceil(vector_length/window_length/overlap) - 2, 2); 

for k=1:window_length
    w(k) = sin((pi * k)/(window_length - 1)); % Sine-window calculation
end



for k=0:ceil(vector_length/overlap/window_length) - 3
    
    index = k * overlap * window_length; % Moving index
    
    X1 = vector(index + 1:index + window_length, 1).*w; 
    X2 = vector(index + 1:index + window_length, 2).*w;
    
    X(:, k+1, 1) = fft(X1); % Output of the first channel
    X(:, k+1, 2) = fft(X2); % Output of the second channel
    
end
end

