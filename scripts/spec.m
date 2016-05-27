window = 1024;
s = spectrogram(drum(:, 1), window, 'yaxis');
hold on
subplot(2,1,1)
plot(s) 
subplot(2,1,2)
plot(drum(:,1));