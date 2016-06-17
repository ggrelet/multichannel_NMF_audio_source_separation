clear all
close all


%% Open sources
[drum, Fd]  = audioread('sources/drum.wav');
[piano, Fp] = audioread('sources/piano.wav');
[voice, Fv] = audioread('sources/voice.wav');
load('sources/mixing_filters_ozerov.mat');
drum=drum(1:2*Fd,1); % 2 premieres secondes
piano=piano(1:2*Fp,1); % idem
voice=voice(1:2*Fv,1); % 2 premieres secondes

%% Convolutive mixing
% from mixin_filter_ozerov.mat
piano_mixed = mix(piano, 1);
drum_mixed = mix(drum, 2);
voice_mixed = mix(voice, 3);

fmix = piano_mixed + drum_mixed + voice_mixed; % convolutive mixture

%% Parameter definition
J=3; % number of instrument
I=2; % stereo
K_partition=[5,5,5];
betaparam=2;
stop=0.0005;
%% Spectrogram
X=spec_cube(fmix,1024,0.5);
X=X(:,1:512,:); % On ne prend pas en compte les hte frequences
F=size(X,2);
N=size(X,3);
s=zeros(J,2*F,N); % 2*F car on va tronquer après
s(1,:,:)=spec(piano,1024,0.5);
s(2,:,:)=spec(drum,1024,0.5);
s(3,:,:)=spec(voice,1024,0.5);
s=s(:,1:512,:);

%% Initialisation de l'algo
W=[];
H=[];
for j=1:J
   [W_temp,H_temp]=nmf_initialization(abs(squeeze(s(j,:,:))),betaparam,stop,K_partition(j)); 
   W=[W W_temp];
   H=[H;H_temp];
end
sigb=zeros(I,I,F);
for f=1:F
    sigb(:,:,f)=0.001*eye(I);    
end


%% Algo sur critÃ¨re d'arret 
[A_est, W_est, H_est, s_est, sigb] = em_step(X, A, W, H, sigb, K_partition);
pause()
% arret = 1;
% while(arret)
%     s_temp = s(1,:,:) + s(2,:,:) + s(3,:,:);
%     [A_est, W_est, H_est, s_est, sigb] = em_step(x, A, W, H, sigb, K_partition);
%     arret = norm(abs(s_step- ( s(1,:,:) + s(2,:,:) + s(3,:,:)))) > 10^-7;
% end
