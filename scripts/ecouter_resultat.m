load('/Users/guillaumegrelet/SIGMA205/functions/signal_estime.mat', 'signal_estime')


A = squeeze(signal_estime(3, :, :));
B = ispec(A, 512);

soundsc(real(B))