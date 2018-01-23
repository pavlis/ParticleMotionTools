function [ sumwtw,sumfw ] = plotslepians( re, im,dt,nstages, f_oversampling )
%plots pairs of slepian wavelets along with time and frequency weight
% functions in 3 subplots.   Last plot produces is overall frequency
% weight in set of nstages divide by 2 decimators. 
% f_oversampling is a multiplier used to zero pad prior to fft to 
% improve appearance of frequency domain plots (has not effect on time
% domain)
% Returns overall time weighting function and frequency weighting
% function
sz=size(re);
nw=sz(1,2);
npts=sz(1,1);
T=npts*dt;
t=-(T-dt)/2:dt:T/2;
df=1/(dt*npts);
df=df/f_oversampling;
fNy=1.0/(2.0*dt);
f= -fNy:df:(fNy-df);
nf=max(size(f));
sumtw=zeros(npts,1);
cw=zeros(nf,1);
cwp=zeros(nf,1);
sumfw=zeros(nf,1);
figure;
for i=1:nw
    % Don't worry about the phase for freq domain result - hence
    % just stick the data in 1 to npts
    cw(1:npts)=complex(re(:,i),im(:,i));
    cwp(1:npts)=complex(im(:,i),re(:,i));
    subplot(3,1,1),plot(t,re(:,i),'b',t,im(:,i),'r');
    xlabel('time (s)');
    title('wavelet');
    subplot(3,1,2),plot(t,abs(cw(1:npts)));
    xlabel('time (s)');
    title('time weighting');
    sumtw=sumtw+abs(cw(1:npts));
    CW=fft(cw);
    CW=fftshift(CW);
    subplot(3,1,3),plot(f,abs(CW));
    xlabel('Frequency (Hz)');
    title('frequency weighting');
    sumfw=sumfw+abs(CW);
    % I think this will remove + and - irregularity (sin-cos term
    % ambiguity)
    CW=fft(cwp);
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
