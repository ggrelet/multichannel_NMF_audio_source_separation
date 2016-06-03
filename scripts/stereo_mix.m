%% Mixture script

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


%% testEm_step script

I=2;
F=1024;
N=1024;
J=4; %nombre de sources
x=zeros(I,F,N); % spectrogramme du signal recu
K=15;
W=zeros(F,K); % NMF estimee de x
H=zeros(K,N); % idem
sigb=zeros(I); %bruit sur les différentes chaines
em_step(x,A,W,H,sigb)