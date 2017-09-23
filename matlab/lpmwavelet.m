function [ w, lambda ] = lpmwavelet( fc, p , dt, M )
%Computed multiwavelet functions, w, using the eigenvector
%equation from Lilly and Park (1995). 
%  fc - center frequency in Hz
%  p - time-bandwidth product (see L&P paper) = fw*M*dt
%  dt - sample interval
%  M - length in samples of wavelets to be generated.  
%
% returns wavelets in Mx4p array w.   lamba is a vector of
% eigenvalues associated with each column of w. 
fw=p/(M*dt)
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
[u,d]=eig(A);
%L&P claim number of useful function is 4*time-bandwidth produce
N=fix(4.0*p);
w=zeros(M,N);
lambda=zeros(N,1);
% eig returns eigenvectors and eigenvalues sorted from smallest to 
% largest - we have reverse that for return 
for i=1:N
    ii=M-i+1;
    w(:,i)=u(:,ii);
    lambda(i)=d(ii,ii)
end

