close all
clear all

% Open sources
[drum, Fd]  = audioread('sources/drum.wav');
[piano, Fp] = audioread('sources/piano.wav');
[voice, Fv] = audioread('sources/voice.wav');
load('sources/mixing_filters_ozerov.mat');

% Get channels
A11 = squeeze(A(1,1,:));
A12 = squeeze(A(1,2,:));
A13 = squeeze(A(1,3,:));

A21 = squeeze(A(2,1,:));
A22 = squeeze(A(2,2,:));
A23 = squeeze(A(2,3,:));

% Make mixed instruments with convolutive coefficients from
% mixin_filter_ozerov.mat
piano_mixed = mix(piano, 1);
drum_mixed = mix(drum, 2);
voice_mixed = mix(voice, 3);

fmix = piano_mixed + drum_mixed + voice_mixed; % Filtered mix
omix = piano + drum + voice; % Original mix
