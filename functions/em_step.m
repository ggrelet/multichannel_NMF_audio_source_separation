function [A_new,W_new,H_new,s,sigb_new] =em_step(x,A,W,H,sigb)

[I,F,N]=size(x); % Nombre de chaines*nombre de fréquence*nombre de tps
[~,J,~]=size(A); % Nombre d'instruments
K_partition=[5 5 5 5];
K_cumsum = cumsum(K_partition);
K=sum(K_partition);
s= %% A INITIALISER
c= %% A INITIALISER

for f=1:F
    Rxx=zeros(I,I);
    Rxs=zeros(I,J);
    Rss=zeros(J,J);
    u=zeros(K,N);
    %% E step
    for n=1:N
        
        
        sigs=zeros(J);
        cpt=1;
        A_ronde=zeros(I,K);
        
        for k=1:K % Calculs concernant les partition de [1,K]
            sigs(cpt,cpt)=sigs(cpt,cpt)+W(f,k)*H(k,n);
            A_ronde(:,k)=A(:,cpt,f);
            if k==K_cumsum(cpt)
                cpt=cpt+1;
            end
        end
        
        sigx=A(:,:,f)*sigs*A(:,:,f)'+sigb(f);
        sigc=diag(W(f,:)*H(:,n));
        Gs=sigs*A(:,:,f)'/(sigx);
        Gc=sigc*A_ronde'/sigx;
        s=Gs*x(:,f,n);
        c=Gc*x(:,f,n);
        Rxx(:,:)=Rxx(:,:)+x(:,f,n)*x(:,f,n)'/N;
        Rxs(:,:)=Rxs(:,:)+x(:,f,n)*s(:,f,n)'/N;
        Rss(:,:)=Rss(:,:)+(s(:,f,n)*s(:,f,n)'+sigs-Gs*A(:,:,f)*sigs)/N;
        u(:,n)=diag(c(:,f,n)*c(:,f,n)'+sigc-Gc*A_ronde*sigc);
        
        
    end
    %% Mstep
    
end

%% M step



end