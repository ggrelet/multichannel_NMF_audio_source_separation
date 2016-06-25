load('signal_estime_1000iter.mat', 'signal_estime')

for k = 1:3
    
    A = squeeze(signal_estime(k , :, :));
    B = ispec(A, 512);

    soundsc(real(B));
    pause(4);
end