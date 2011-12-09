function MRSplot(freq, gabadata)
%function MRSplot(freq, gabadata)
% Plots MRS data loaded by LoadPfileGABA

%figure

plot(freq, real(gabadata), 'k');
set(gca,'XDir','reverse');
oldaxis = axis;
axis([0 5 0.5*oldaxis(3) oldaxis(4)])
