nbins = 500;
vect_bins= linspace(min(Xfault),max(Xfault),nbins);
Y_bins = zeros(1,nbins);

for i = 2:1:nbins-1
    Y_bins(i) = mean(Yfault(Xfault> vect_bins(i-1) & Xfault< vect_bins(i+1)));
end

vect_bins = vect_bins(Y_bins>0);
Y_bins = Y_bins(Y_bins>0);

figure
plot(Xfault,Yfault,'k.')
hold on
plot(vect_bins(Y_bins>0),Y_bins(Y_bins>0),'r--')
