I=2;
F=1024;
N=65536;
J=4; %nombre de sources
x=zeros(I,F,N); % spectrogramme du signal reçu
A=zeros(I,J,F); % Filtre qui fait le mixage
K=15;
W=zeros(F,K); % NMF estimee de x
H=zeros(K,N); % idem
s=zeros(J,F,N); % spectrogramme du signal des sources
sigb=zeros(I); %bruit sur les différentes chaines
em_step(x,A,W,H,s,sigb)
