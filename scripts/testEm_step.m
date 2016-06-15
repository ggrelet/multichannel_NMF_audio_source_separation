%% SCRIPT DE TEST DE LA FONCTION EM_STEP
% Cette fonction teste l'algorithme EM avec des valeurs simples et courtes.
% Son temps de réponse est de quelques secondes là où il aurait fallu
% plusieures minutes sur les signaux drum, piano et voice.
% La première étape fixe les paramètres tandis que la seconde lance
% l'algorithme à proprement parler.
%% Parametres

I=2; % Nombre de chaines
F=16;
N=50; % Nombre de points et donc d'itérations de l'algorithmes
J=3; % Nombre de sources
K_partition=[5 4 5];
K=sum(K_partition);
x=randn(I,F,N); % Spectrogramme du signal
A=randn(I,J,F); % Filtre qui fait le mixage

W=rand(F,K); % Les matrices doivent être positives pour que l'algorithme fonctionne
H=rand(K,N); 

sigb=zeros(I,I,F);
for f=1:F
    sigb(:,:,f)=0.001*eye(I);    
end

re=zeros(1,N);
im=re;

%% Test
for i=1:N
   [A_new,W_new,H_new,~,sigb_new,criterion]=em_step(x,A,W,H,sigb,K_partition);
   A=A_new;
   W=W_new;
   H=H_new;
   re(i)=real(criterion);
   im(i)=imag(criterion);
   sigb=sigb_new;
end

%% Résultats
% Les résultats suivants nous montrent plusieures choses.
% La première est que même si elle est très petite, le critère possède une
% partie imaginaire non nulle. La seconde est que la partie réelle est bien décroissante. Enfin,
% le résultat affiché res nous indique que cette décroissance se fait de
% manière strictement monotone.
close all
figure
hold on
subplot(2,1,1)
plot(re)
title('Partie réelle')
subplot(2,1,2)
plot(im, 'r')
title('Partie imaginaire')
hold off

d=diff(re);
res = sum(d>0);
display(res);
