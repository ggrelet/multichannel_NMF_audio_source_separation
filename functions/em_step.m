function [A_new,W_new,H_new,s,sigb_new] =em_step(x,A,W,H,sigb)


%%
% Une interation de EM
[I,F,N]=size(x); % Nombre de chaines*nombre de fréquence*nombre de tps
[~,J,~]=size(A); % Nombre d'instruments
K_partition=[5 5 5];
K_cumsum = cumsum(K_partition);
K=sum(K_partition);
s= zeros(J,F,N);
c= zeros(K,F,N); %% Pas sur
A_new = zeros(size(A));
sigb_new = zeros(size(sigb));
u=zeros(K,F,N);
H_new=zeros(size(H));
W_new=zeros(size(W));

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
        
        sigx=A(:,:,f)*sigs*A(:,:,f)'+sigb(:,:,f);
        sigc=diag(diag(W(f,:)*H(:,n)));
        Gs=sigs*A(:,:,f)'/(sigx);
        Gc=sigc*A_ronde'/sigx;
        s(:,f,n)=Gs*x(:,f,n);
        c(:,f,n)=Gc*x(:,f,n);
        Rxx(:,:)=Rxx(:,:)+x(:,f,n)*x(:,f,n)'/N;
        Rxs(:,:)=Rxs(:,:)+x(:,f,n)*s(:,f,n)'/N;
        Rss(:,:)=Rss(:,:)+(s(:,f,n)*s(:,f,n)'+sigs-Gs*A(:,:,f)*sigs)/N;
        u(:,f,n)=diag(c(:,f,n)*c(:,f,n)'+sigc-Gc*A_ronde*sigc);
        
    end
    %% Mstep
    temp=Rxs/Rss;
    A_new(:,:,f)=temp;
    
    sigb_new(:,:,f)=diag(diag(Rxx-temp*Rxs'-Rxs*temp'+temp*Rss*temp'));
    for k=1:K
        W_new(f,k)=sum(squeeze(u(k,f,:)).'./H(k,:))/N;
    end
    
end
for n=1:N
    for k=1:K
        H_new(k,n)=sum(squeeze(u(k,:,n)).'./W_new(:,k))/F;
    end
end

%% Normalization of A_new W_new and H_new
D=zeros(J,J,F);
for f=1:F
    
    for j=1:J
        D(j,j,f)=sqrt(sum(abs(A_new(:,j,f).^2)))*exp(1i*angle(A_new(1,j,f)));
    end
    A_new(:,:,f)=A_new(:,:,f)/D(:,:,f);
    % on a ainsi la sommes sur i  A_ij,f = 1 et A_1j reel > 0
end
ind=[1 K_cumsum]; %indices des H_j et K_j
for j=1:J
    cardKj=ind(j+1)-ind(j);
    lambda=zeros(cardKj,cardKj);
    % Attention W deux dim pas une
    W(:,ind(j):ind(j+1))=diag(abs(squeeze(D(j,j,:))).^2)*W(:,ind(j):ind(j+1));
    for k=ind(j):ind(j+1)
        lambda(k-ind(j)+1,k-ind(j)+1)=sum(W(:,k));
    end
    W(:,ind(j):ind(j+1))=W(:,ind(j):ind(j+1))*lambda;
    H(ind(j):ind(j+1),:)=lambda*H(ind(j):ind(j+1),:);
    % multiplier J_j par lambda
end

end