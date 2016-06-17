function [A_new,W_new,H_new,s,sigb_new,criterion] = em_step(x,A,W,H,sigb,K_partition)


%%
% Une interation de EM
[I,F,N]=size(x); % Nombre de chaines*nombre de fréquence*nombre de tps
[~,J,~]=size(A); % Nombre d'instruments
if J~=size(K_partition,2)
    error('Attention K_partition n''a pas J elements')
end
K_cumsum = cumsum(K_partition);
ind=[0 K_cumsum]; %indices des H_j et K_j
K=sum(K_partition);
matrice=zeros(J,K); % matrice pour calculer A_ronde
for j=1:J
   matrice(j,ind(j)+1:ind(j+1))=ones(1,K_partition(j)); 
end

s= zeros(J,F,N);
c= zeros(K,F,N);
A_new = zeros(size(A));
sigb_new = zeros(size(sigb));
u=zeros(K,F,N);
H_new=zeros(size(H));
W_new=zeros(size(W));
sigx=zeros(I,I,F,N);

for f=1:F
   %% E Step
   % Rxx=zeros(I,I);
    Rss=zeros(J,J);
    A_ronde=A(:,:,f)*matrice;   
    sigs_n=zeros(J,N); % Matrice diagonale, on ne stocke que la diagonale
    for j=1:J
        sigs_n(j,:)=W(f,ind(j)+1:ind(j+1)) * H(ind(j)+1:ind(j+1),:);
    end
    for n=1:N
        sigx(:,:,f,n)=A(:,:,f)*diag(sigs_n(:,n))*A(:,:,f)'+sigb(:,:,f);
        sigc=(diag(W(f,:).*H(:,n).'));
        Gs=diag(sigs_n(:,n))*A(:,:,f)'/(sigx(:,:,f,n)); % calculer avant l'inverse
        Gc=sigc*A_ronde'/sigx(:,:,f,n);
        s(:,f,n)=Gs*x(:,f,n);
        c(:,f,n)=Gc*x(:,f,n); 
        Rss=Rss+(s(:,f,n)*s(:,f,n)'+diag(sigs_n(:,n))-Gs*A(:,:,f)*diag(sigs_n(:,n)))/N;
        % u(:,:,f)=abs(c(:,:,f)).^2 + sigc(:,:,f) ...
        % - permute(sum(bsxfun(@times, Gc(:,:,:,f), A_ronde(:,:,f).'),2),[1,3,2]) .* sigc(:,:,f);
        u(:,f,n)=diag(c(:,f,n)*c(:,f,n)'+sigc-Gc*A_ronde*sigc);

    end
    
    x=permute(x,[1 3 2]);
    s=permute(s,[1 3 2]);
    Rxx=x(:,:,f)*x(:,:,f)'/N;
    Rxs=x(:,:,f)*s(:,:,f)'/N;
    x=permute(x,[1 3 2]);
    s=permute(s,[1 3 2]);
    %% Mstep
    temp=Rxs/Rss;
    A_new(:,:,f)=temp;
    
    sigb_new(:,:,f)=diag(diag(Rxx-temp*Rxs'-Rxs*temp'+temp*Rss*temp'));
    %u2 = permute(u,[1,3,2]);
    %u_temp = u2(:,:,f);
    u_temp=permute(u(:,f,:),[1,3,2]);
    for k=1:K
        W_new(f,k)=sum(u_temp(k,:)./H(k,:))/N;
    end
    
end
for k=1:K
    u_temp=squeeze(u(k,:,:));
    for n=1:N
        H_new(k,n)=sum(u_temp(:,n)./W_new(:,k))/F;
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

for j=1:J
    
    lambda=zeros(K_partition(j));
    W(:,ind(j)+1:ind(j+1))=diag(abs(squeeze(D(j,j,:))).^2)*W(:,ind(j)+1:ind(j+1));
    for k=ind(j)+1:ind(j+1)
        lambda(k-ind(j),k-ind(j))=sum(W(:,k));
    end
    W(:,ind(j)+1:ind(j+1))=W(:,ind(j)+1:ind(j+1))/lambda;
    H(ind(j)+1:ind(j+1),:)=lambda*H(ind(j)+1:ind(j+1),:);
end

%% Calcul du critere v
v=0;
for n=1:N
    for f=1:F
        v=v+trace((x(:,f,n)*x(:,f,n)')/sigx(:,:,f,n))+log((det(sigx(:,:,f,n))));
    end
end
criterion=v;
end