function MRSplotspec(MRS_struct)
%function MRSplotspec(MRS_struct)
% Plots MRS data loaded by LoadPfileGABA

%figure

plot(MRS_struct.freq, real(MRS_struct.FullDataOFF));
set(gca,'XDir','reverse');
oldaxis = axis;
axis([0 4.5  0.1*oldaxis(3) 0.25*oldaxis(4)])
