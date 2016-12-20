function [vector] = ispec(X, window_length)
%ispec fonction spectrogramme inverse

	vector = [];

	w = zeros(window_length, 1);
	for k=1:window_length
	    w(k) = sin((pi * k)/(window_length - 1)); % Sine-window calculation
	end

	for k = 1:2:size(X, 2)
	  
	    A = ifft(X(:, k));
	    A = A./w;
	    vector = [vector ; A];
	end
end