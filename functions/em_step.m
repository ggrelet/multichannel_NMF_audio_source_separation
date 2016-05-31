function [A_new,W_new,H_new,s_new,sigb_new] =em_step(x,A,W,H,s,sigb)

    [I,F,N]=size(x); % Nombre de chaines*nombre de fréquence*nombre de tps
    [~,J,~]=size(A); % Nombre d'instruments
    K_partition=[5 5 5 5];
    K_cumsum = cumsum(K_partition);
    K=sum(K_partition);
    
    %% E step
    Rxx=zeros(I,I,F);
    Rxs=zeros(I,J,F);
    Rss=zeros(J,J,F);
    u=zeros(K,F,N)
    for n=1:N
        for f=1:F
            sigs=zeros(J);
            sigx=A(:,:,f)*sigs*A(:,:,f)'+sigb(f);
            sigc=diag(W(f,:)*H(:,n));
            Gs=sigs*A(:,:,f)'/(sigx);
            Rxx(:,:,f)=Rxx(:,:,f)+x(:,f,n)*x(:,f,n)';
            Rxs(:,:,f)=Rxs(:,:,f)+x(:,f,n)*s(:,f,n)';
            Rss(:,:,f)=Rss(:,:,f)+s(:,f,n)*s(:,f,n)'+sigs-Gs*A(:,:,f)*sigs;
            
            sigc=...
            A_ronde=...
            u(:,f,n)=diag(c(:,f,n)*c(:,f,n)'+sigc-Gc*A_ronde*sigc);
            
            
            sigs=zeros(J);
            cpt=1;
            for k=1:K
                sigs(cpt,cpt)=sigs(cpt,cpt)+W(f,k)*H(k,n);
                if k==K_cumsum(cpt)
                    cpt=cpt+1;
                end
            end
        end
    end
    Rxx=Rxx/N;
    Rxs=Rxs/N;
    Rss=Rss/N;
    
    %% M step
    


end