clear all
close all


%% Open sources
[drum, Fd]  = audioread('sources/drum.wav');
[piano, Fp] = audioread('sources/piano.wav');
[voice, Fv] = audioread('sources/voice.wav');
load('sources/mixing_filters_ozerov.mat');
drum=drum(:,1);
piano=piano(:,1);
voice=voice(:,1);

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
K_param=[5,5,5];
betaparam=2;
stop=10^-7;
%% Initialisation de l'algo
[W,H]=nmf_initialization(piano,betaparam,stop);
[W_temp,H_temp]=nmf_initialization(drum,betaparam,stop);
W=[W W_temp];
H=[H;H_temp];
[W_temp,H_temp]=nmf_initialization(voice,betaparam,stop);
W=[W W_temp];
H=[H;H_temp];


%% Algo sur critère d'arret 
[A, W, H, s, sigb] = em_step(x, A, W, H, sigb, K_partition);
arret = 1;
while(arret)
    s_temp = s(1,:,:) + s(2,:,:) + s(3,:,:);
    [A, W, H, s, sigb] = em_step(x, A, W, H, sigb, K_partition);
    arret = norm(abs(s_step- ( s(1,:,:) + s(2,:,:) + s(3,:,:)))) > 10^-7;
end
