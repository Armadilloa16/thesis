function [Shat,nCal,meanCalIntensity,CV] = normaliseDataset(L,LXY,Vars,calibrantMZ,tol,var,includeCV,ppmTol)
% Normalisation to internal calibrants

[isOctave,matType,matExt] = checkIsOctave();

if nargin < 5
    tol = 50;
end
if nargin < 6
   var = 'intensity';
end
if nargin < 7
    includeCV = true(1,1);
end
if nargin < 8
    ppmTol = true;
end

CV = struct();
CV.calibrantMZ = calibrantMZ;
CV.tol = tol;
CV.ppmTol = ppmTol;

% Initialise matrix for storing calibrant intensities.
calibrantIntensities = zeros(length(calibrantMZ),length(L));
    
% Similar to the sortPeaks function.
% Compile all peaklists into a comprehensive list for easy access.
includeVar = strcmp(Vars,var);
allPeaks = zeros(sum(LXY(:,3)),4);
curLoc = 1;
for spec = 1:size(LXY,1)
    nPeaks = LXY(spec,3);
    allPeaks(curLoc:(curLoc+nPeaks-1),1) = spec*ones(LXY(spec,3),1);
    allPeaks(curLoc:(curLoc+nPeaks-1),2) = 1:nPeaks;
    allPeaks(curLoc:(curLoc+nPeaks-1),3) = L{spec}(:,1);
    % Here the difference to the sortPeaks function -- include
    % intensities.
    allPeaks(curLoc:(curLoc+LXY(spec,3)-1),4) = L{spec}(:,includeVar);
    curLoc = curLoc + nPeaks;
end

% Now, for each calibrant, we will extract intensities.
for cal_idx = 1:length(calibrantMZ)
    % Look within a plus/minus tol Da window of the true calibrant
    % mass for peaks.
    if ppmTol
        subset_bin = abs(allPeaks(:,3)-calibrantMZ(cal_idx)) <= tol*calibrantMZ(cal_idx)/1000000;
    else
        subset_bin = abs(allPeaks(:,3)-calibrantMZ(cal_idx)) <= tol;
    end
    unique_spectra = unique(allPeaks(subset_bin,1));
    % Check to see if there are any spectra which contain more than 1
    % peak within that window.
    if length(unique_spectra) ~= sum(subset_bin)
        % and if there are, choose the peak with m/z closest to the
        % true value to represent the calibrant, and print a warning.
        for i = unique_spectra(histc(allPeaks(subset_bin,1),unique_spectra) > 1)
            dup_peaks_logic = (subset_bin & allPeaks(:,1) == i);
            [~,I] = min(abs(allPeaks(dup_peaks_logic,3) - calibrantMZ(cal_idx)));
            subset_bin(dup_peaks_logic) = (1:sum(dup_peaks_logic) == I);
        end
        % I deal with this by editing subset_bin here to only include
        % the closest peak in each spectrum.
        disp('warning: duplicate peaks present -- closest peak (in m/z) chosen.')
    end
    unique_intensitites = zeros(length(L),1);
    unique_intensitites(allPeaks(subset_bin,1)) = allPeaks(subset_bin,4);
    % Store the results in calibrantIntensities. This will be our data
    % for the normalisation regression.
    calibrantIntensities(cal_idx,:) = unique_intensitites;
end
    
% Regression 
X = zeros(size(calibrantIntensities));
X(calibrantIntensities>0) = log(calibrantIntensities(calibrantIntensities>0));
nCal = sum(X>0);
Xbar = zeros(size(X,1),1);
for i = 1:size(X,1)
    Xbar(i) = mean(X(i,X(i,:) > 0));
end
Shat = zeros(1,size(X,2));
for j = find(nCal > 0)
    Shat(j) = mean(X(X(:,j)>0,j) - Xbar(X(:,j)>0));
end
meanCalIntensity = mean(Xbar);

% CV
if includeCV
    Shat_cv = cell(1,size(X,1)-1);

    choiceL = cell(1,length(Shat_cv));
    for n_calibrants = 1:length(Shat_cv)
        temp = perms([ones(1,n_calibrants) zeros(1,size(calibrantIntensities,1)-n_calibrants)]);
        choiceL{n_calibrants} = cell(1,size(temp,1));
        for i = 1:size(temp,1)
            choiceL{n_calibrants}{i} = num2str(temp(i,:));
        end
        choiceL{n_calibrants} = unique(choiceL{n_calibrants});
        for i = 1:length(choiceL{n_calibrants})
            choiceL{n_calibrants}{i} = str2num(choiceL{n_calibrants}{i});
        end
    end

    for n_calibrants = 1:length(Shat_cv)
        Shat_cv{n_calibrants} = zeros(nchoosek(size(calibrantIntensities,1),n_calibrants),size(calibrantIntensities,2));
        for i = 1:length(choiceL{n_calibrants})
            temp = calibrantIntensities(choiceL{n_calibrants}{i}>0,:);
            X = temp;
            X(temp>0) = log(temp(temp>0));
            Xbar = zeros(size(X,1),1);
            for k = 1:size(X,1)
                Xbar(k) = mean(X(k,X(k,:) > 0));
            end

            for j = 1:size(X,2)
                if sum(X(:,j)>0) > 0
                    Shat_cv{n_calibrants}(i,j) = mean(X(X(:,j)>0,j) - Xbar(X(:,j)>0));
                end
            end
        end
    end
    CV.Shat = Shat_cv;
    CV.choiceL = choiceL;
end

end


