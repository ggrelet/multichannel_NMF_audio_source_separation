load('signal_estime_100iter_1024freq.mat', 'signal_estime')


A = squeeze(signal_estime(1, :, :));
B = ispec(A, 1024);

soundsc(real(B),16000)