clear all
I=2;
F=1024;
N=1024;
J=3; %nombre de sources
x=zeros(I,F,N); % spectrogramme du signal recu
A=zeros(I,J,F); % Filtre qui fait le mixage
K=15;
W=zeros(F,K); % NMF estimee de x
H=zeros(K,N); % idem
sigb=zeros(I,I,F); %bruit sur les différentes chaines
em_step(x,A,W,H,sigb);
