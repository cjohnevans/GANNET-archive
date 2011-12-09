function MRS_struct_new = MRSremovespectrum(MRS_struct,pptNum)
% function MRS_struct_new = MRSremovespectrum(MRS_struct,pptNum) 
% remove spectrum from MRS structure
% if removing multiple spectra from a dataset, start with the last in the 
% struct to avoid confusion!


rmv=input(['    Remove  ' MRS_struct.pfile{pptNum} '? (y/n) ' ], 's');

if strcmp(rmv,'y')
    disp( ['    Removing ' MRS_struct.pfile{pptNum}]);
    MRS_struct_new = MRS_struct;
    MRS_struct_new.pfile(pptNum) = [];
    MRS_struct_new.Navg(pptNum) = [];
    MRS_struct_new.waterspec(pptNum,:) = [];
    MRS_struct_new.waterfreq(pptNum,:) = [];
    MRS_struct_new.fwhmHz(pptNum) = [];
    MRS_struct_new.gabanoalign(pptNum,:) = [];
    MRS_struct_new.gabaspec(pptNum,:) = [];
    MRS_struct_new.sumspec(pptNum,:) = [];
    MRS_struct_new.FreqStdevHz(pptNum) = [];
    MRS_struct_new.phase(pptNum) = [];
    MRS_struct_new.phase_firstorder(pptNum) = [];
    MRS_struct_new.Rejects(pptNum) = [];
    MRS_struct_new.GABA_SNR(pptNum) = [];
    MRS_struct_new.gabaArea(pptNum) = [];
    MRS_struct_new.GABAFWHM(pptNum) = [];
    MRS_struct_new.LGModelParams(pptNum,:) = [];
    MRS_struct_new.waterArea(pptNum) = [];
    MRS_struct_new.GABAModelFit(pptNum) = [];
    MRS_struct_new.gabaspec_scaled(pptNum,:) = [];
    MRS_struct_new.gabaiu(pptNum) = [];
else
    disp( ['    Doing nothing.']);
    MRS_struct_new = -1;
end

