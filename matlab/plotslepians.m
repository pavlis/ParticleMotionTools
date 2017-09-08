function  plotslepians( re, im )
%plots pairs of slepian wavelets with on plot per pair
sz=size(re);
nw=sz(1,2);
npts=sz(1,1);
for i=1:nw
    plot(re(:,i));
    hold;
    plot(im(:,i));
    %if( i~= nw) 
     %   figure;
    %end
    figure;
end
for i=1:nw
    S=fft(re(:,i));
    semilogy(abs(S(1:npts/2)));
    if (i == 1) 
        hold;
    end
end
figure;
Ssum=zeros(npts,1);
for i=1:nw
    S=fft(im(:,i));
    semilogy(abs(S(1:npts/2)));
    if (i == 1)
        hold;
    end
    Ssum=Ssum+abs(S);
end
figure;
semilogy(Ssum);
