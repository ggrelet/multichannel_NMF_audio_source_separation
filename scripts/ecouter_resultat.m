<<<<<<< HEAD
load('signal_estime_100iter_1024freq.mat', 'signal_estime')
=======
load('signal_estime_1000iter.mat', 'signal_estime')
>>>>>>> 797e7041b0e1b07632b0450f000f9f4d169491f7

for k = 1:3
    
    A = squeeze(signal_estime(k , :, :));
    B = ispec(A, 512);

<<<<<<< HEAD
A = squeeze(signal_estime(1, :, :));
B = ispec(A, 1024);

soundsc(real(B),16000)
=======
    soundsc(real(B));
    pause(4);
end
>>>>>>> 797e7041b0e1b07632b0450f000f9f4d169491f7
