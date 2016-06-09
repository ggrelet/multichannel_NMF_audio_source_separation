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

fmix = piano_mixed + drum_mixed + voice_mixed; % Filtered mix
omix = piano + drum + voice; % Instantaneous mix

%% Spectrogram
X=spec(fmix,1024,0.5);

%% Initialisation de l'algo

%% Algo sur critère d'arret 

