function [A_new,W_new,H_new,s_new] =em_step(x,A,W,H,s)

    [I,F,N]=size(x);
    Rxx=zeros(I,I,F);
    Rxs=zeros()
    for n=1:N
        for f=1:F
            Rxx(:,:,f)=Rxx(:,:,f)+x(:,f,n)*x(:,f,n)';
        end
    end
    Rxx=Rxx/N;
    
    


end