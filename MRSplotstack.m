function MRSplotstack(MRS_struct)
%function MRSplotstack(MRS_struct)
% Plots MRS data loaded by MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA

numspec = length(MRS_struct.gabaspec(:,1));

% Find Water amplitude max, across all Pfiles
waterheight = abs(max(MRS_struct.waterspec,[],2));
heightrescale = repmat((1./waterheight), [1 length(MRS_struct.gabaspec(1,:))]);
SpectraToPlot = MRS_struct.gabaspec .* heightrescale;

% Estimate baseline from between Glx and GABA
specbaseline = mean(real(SpectraToPlot(:,17250:17650)),2);

% averaged gaba height across all scans - to estimate stack spacing
gabaheight = abs(max(SpectraToPlot(:,17250:18000),[],2));
gabaheight = mean(gabaheight);

plotstackoffset = [ 0 : (numspec-1) ]';
plotstackoffset = plotstackoffset * gabaheight;
plotstackoffset = plotstackoffset - specbaseline;

SpectraToPlot = SpectraToPlot + ...
    repmat(plotstackoffset, [ 1  length(MRS_struct.gabaspec(1,:))]);

%figure(99)
plot(MRS_struct.freq, real(SpectraToPlot));
legendtxt = regexprep(MRS_struct.pfile, '_','-');

legend(legendtxt)
set(gca,'XDir','reverse');
oldaxis = axis;
% yaxis max = top spec baseline + 2*meangabaheight
yaxismax = (numspec + 2) *gabaheight; % top spec + 2* height of gaba
yaxismin =  - 2* gabaheight; % extend 2* gaba heights below zero

axis([0 5  yaxismin yaxismax])

