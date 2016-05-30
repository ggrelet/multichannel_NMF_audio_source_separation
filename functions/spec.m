function [X] = spec(vector, window_length, overlap)
%SPEC Spectrogram function with sine-window
%   vector : input x vector
%   window_length : size of the sine-window, typically 1024 samples
%   overlap : overlapping coefficient, typically 50%

vector_length = size(vector, 1); % Size of the input vector (the sound)

w = zeros(window_length, 1); % Window vector
X = zeros(window_length, ceil(vector_length/window_length/overlap)); 

for k=1:window_length
    w(k) = sin((pi * k)/(window_length - 1)); % Sine-window
end



for k=0:ceil(vector_length/overlap/window_length) 
    index = k * overlap * window_length; % Change index
    X1 = vector(index + 1:index + window_length).*w;
    X(:, k+1) = fft(X1); % Output
end
end

