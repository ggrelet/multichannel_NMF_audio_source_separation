clear all

%% Parametres

I=2; % nombre de chaines
F=16;
N=20;
J=3; %nombre de sources
K=15;
K_partition=[5 5 5];
x=randn(I,F,N); % spectrogramme du signal
A=randn(I,J,F); % Filtre qui fait le mixage

W=rand(F,K);
H=rand(K,N);
sigb=zeros(I,I,F);
for f=1:F
    sigb(:,:,f)=0.001*eye(I);    
end

%% Test
for i=1:5
   [A_new,W_new,H_new,~,sigb_new]=em_step(x,A,W,H,sigb,K_partition);
   A=A_new;
   W=W_new;
   H=H_new;
   sigb=sigb_new;
end
