function [A_new,W_new,H_new,s_new,sigb_new] =em_step(x,A,W,H,s,sigb)

    [I,F,N]=size(x); % Nombre de chaines*nombre de fréquence*nombre de tps
    [~,J,~]=size(A); % Nombre d'instruments
    K=15; %ordre de la NMF
    Rxx=zeros(I,I,F);
    Rxs=zeros(I,J,F);
    Rss=zeros(J,J,F);
    for n=1:N
        for f=1:F
            sigs=zeros(J);
            sigx=A(:,:,f)*sigs*A(:,:,f)'+sigb(f);
            sigc=diag(W(f,:)*H(:,n));
            Gs=sigs*A(:,:,f)/(sigx);
            Rxx(:,:,f)=Rxx(:,:,f)+x(:,f,n)*x(:,f,n)';
            Rxs(:,:,f)=Rxs(:,:,f)+x(:,f,n)*s(:,f,n)';
            Rss(:,:,f)=Rss(:,:,f)+s(:,f,n)*s(:,f,n)'+sigs-Gs*A(:,:,f)*sigs;
            sigs(:,:,f,n)=diag(sum(W(f,:).*H(:,n)));
        end
    end
    Rxx=Rxx/N;
    Rxs=Rxs/N;
    Rss=Rss/N;
    


end