close all
clear all

[drum, Fd]  = audioread('sources/drum.wav');
[piano, Fp] = audioread('sources/piano.wav');
[voice, Fv] = audioread('sources/voice.wav');
load('sources/mixing_filters_ozerov.mat');

A11 = squeeze(A(1,1,:));
A12 = squeeze(A(1,2,:));
A13 = squeeze(A(1,3,:));

A21 = squeeze(A(2,1,:));
A22 = squeeze(A(2,2,:));
A23 = squeeze(A(2,3,:));

piano_mixed = mix(piano, 1);
drum_mixed = mix(drum, 2);
voice_mixed = mix(voice, 3);

fmix = piano_mixed + drum_mixed + voice_mixed;
omix = piano + drum + voice;
