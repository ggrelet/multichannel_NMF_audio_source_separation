clear all
close all


%% Open sources
[drum, Fd]  = audioread('sources/drum.wav');
[piano, Fp] = audioread('sources/piano.wav');
[voice, Fv] = audioread('sources/voice.wav');
load('sources/mixing_filters_ozerov.mat');


%% Convolutive mixing
% from mixin_filter_ozerov.mat
piano_mixed = mix(piano, 1);
drum_mixed = mix(drum, 2);
voice_mixed = mix(voice, 3);

fmix = piano_mixed + drum_mixed + voice_mixed; % convolutive mixture

%% Spectrogram
X=spec(fmix,1024,0.5);

%% Parameter definition
J=3; % number of instrument
I=2; % stereo
%% Initialisation de l'algo

for j=1:J
    [W,H]=nmf_initialization()
end
%% Algo sur critère d'arret 
 em_step(x,A,W,H,sigb)
while(arret)
    
    em_step(x,A,W,H,sigb)
end
