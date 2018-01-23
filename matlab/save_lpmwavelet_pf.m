function [ sumtw,sumfw ] = save_lpmwavelet_pf( u, fc, p, dt, nstages, f_oversampling, fname )
%Creates a parameter file from multiwavlets stored in u.   It is
% assumed the wavelets are adjacentpairs in columns of u that can
% be placed in real and imaginary parts.   Always plots u as 
% a quality control.   
% Other inputs:
% fc - center frequency used for generating wavelets
% p - time-bandwidth product of these wavelets
% dt - sample interval used with fc in original calculation of wavelets
% (note - the point is fc, p, and dt should be the same as those sent
% to lpmwavelet.m script.)
% nstages - has no effect on pf output, but displays series of frequency
%   domain plots for nstages of divide by 2 decimation
% f_oversampling - controls zero padding for frequency domain plots.  
%   also has not effect on pf output, but improves the appearance of
%   frequency domain plots
% fname - output file name (string)
%plotlpmwavelet(u);
%title('multiwavelet pairs being saved');
%figure;
[m,n]=size(u);
T=m*dt;
t=-(T-dt)/2:dt:T/2;
nw=n/2;
for i=1:nw
    ii=2*(i-1)+1;
    cw(:,i)=complex(u(:,ii),u(:,ii+1));
    cwp(:,i)=complex(u(:,ii+1),u(:,ii));
end
fNy=1.0/(2.0*dt);
df=1/T;
df=df/f_oversampling;
fNy=1.0/(2.0*dt);
f= -fNy:df:(fNy-df);
nf=max(size(f));
sumtw=zeros(m,1);
sumfw=zeros(nf,1);
cws=zeros(nf:1);

f=-fNy:df:(fNy-df);
for i=1:nw
    
    subplot(3,1,1),plot(t,real(cw(:,i)),'b',t,imag(cw(:,i)),'r');  
    xlabel('time (s)');
    title('wavelet');
    subplot(3,1,2),plot(t,abs(cw(:,i)));
    xlabel('time (s)');
    title('time weighting');
    sumtw=sumtw+abs(cw(:,1));
    cws=zeros(nf,1);
    cws(1:m)=cw(:,i);
    CW=fft(cws);
    CW=fftshift(CW);
    subplot(3,1,3),plot(f,abs(CW));
    xlabel('Frequency (Hz)');
    title('frequency weighting');
    sumfw=sumfw+abs(CW);
    % I think this will remove + and - irregularity (sin-cos term
    % ambiguity)
    cws=zeros(nf,1);
    cws(1:m)=cw(:,i);
    CW=fft(cws);
    CW=fftshift(CW);
    sumfw=sumfw+abs(CW);
    figure;
end
plot(f,sumfw);
hold;
for i=1:nstages-1
    f=f/2;
    plot(f,sumfw);
end
xlabel('Frequency (Hz)');
title('Overall frequency weighting');
fid=fopen(fname,'w');
fprintf(fid,'nsamples    %d\n',m);
fprintf(fid,'nwavelets   %d\n',nw);
fNy=1.0/(2*dt);
f0=fc/fNy;
fw=p/(m*dt);
fwnd=fw/fNy;
fprintf(fid,'f0  %f\n',f0);
fprintf(fid,'fw  %f\n',fwnd);
fprintf(fid,'wavelets &Tbl{\n');
re=real(cw);
im=imag(cw);
for j=1:nw
	for i=1:m
		fprintf(fid,'%20.16f %20.16f\n',re(i,j),im(i,j));
	end
end
fprintf(fid,'}\n');
fclose(fid);
