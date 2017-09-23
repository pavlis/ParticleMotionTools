function [ w, lambda ] = mwavlet( fc, fw, dt, M )
%Computed multiwavelet functions, w, using the eigenvector
%equation from Lilly and Park (1995).
A=zeros(M,M);
for n=1:M
    for m=1:M
        a0=2.0*pi*dt*(m-n);
        if(n==m)
            %This is an assymptotic form easily derived by L'Hospital's
            %rule not mentioned in Lilly and Park or Lorie Bear's
            %dissertation.  dt factor took me a bit to work out
            A(m,n)=4.0*fw*dt;
        else
            A(m,n)=(sin(a0*(fc+fw))-sin(a0*(fc-fw)))/(pi*(m-n));
        end
    end
end
[w,lambda]=eig(A);
