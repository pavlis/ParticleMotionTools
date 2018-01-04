function  plotslepians( re, im,tstr )
%plots pairs of slepian wavelets with one plot per pair
sz=size(re);
nw=sz(1,2);
npts=sz(1,1);
df=1/npts;
f=0.0:df:0.5;   % to plot to Nyquist
nf=max(size(f));
figure;
for i=1:nw
    plot(re(:,i));
    hold;
    plot(im(:,i));
    %if( i~= nw) 
     %   figure;
    %end
    title(tstr);
    figure;
end
for i=1:nw
    S=fft(re(:,i));
    semilogy(f,abs(S(1:nf)));
    if (i == 1) 
        hold;
    end
end
title(tstr);
figure;
Ssum=zeros(npts,1);
for i=1:nw
    S=fft(im(:,i));
    semilogy(f,abs(S(1:nf)));
    if (i == 1)
        hold;
    end
    Ssum=Ssum+abs(S);
end
title(tstr);
figure;
semilogy(f,Ssum(1:nf));
title(tstr);
