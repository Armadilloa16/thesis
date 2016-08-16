

% Read peaklists.
readPeaklistsDriver_Vtma

% Binning
binPeaksDriver_Vtma

% Binary Spatial Smooth
spatialSmoothDriver_Vtma

% Produce binned non-binary data
binNonBinaryDataDriver_Vtma

% Normalisation
normalisationDriver_Vtma

% Annotation
manualClassAnnotations_Vtma

% Summarise Regions
regionSummaryDriver_Vtma

% Summarise Patients
patientSummaryDriver_Vtma

% PCA on patient data
ppcaDriver_Vtma

% Array of classification attempts
summaryClassificationDriver_Vtma

% Print output for use in knitr
% Majority (wiggle) classification 
majorityClassificationDriver_Vtma
% Details on which patients were 
%   classified what in each wiggle.
classificationDetailsDriver_Vtma
% Variables selected in the cca
%   variable selection step.
ccaVariablesDriver_Vtma