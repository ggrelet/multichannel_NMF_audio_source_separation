function [A_new,W_new,H_new,s_new] =em_step(x,A,W,H,s,sigb)

    [I,F,N]=size(x);
    Rxx=zeros(I,I,F);
    Rxs=zeros(I,J,F);
    Rss=zeros(J,J,F);
    sigs=zeros(J,J,F,N);
    Gs=zeros(J,J,F,N);
    for n=1:N
        for f=1:F
            sigs=
            sigx=A(:,:,f)*sigs*A(:,:,f)'+sigb(f)
            Gs=sigs*A(:,:,f)*inv(sigx)
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