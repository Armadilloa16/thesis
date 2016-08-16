function [mbdata,vbincentrs,fExists,emptySpec,nIter,converged] = spatiallySmooth(mbdata,vbincentrs,fExists,emptySpec,smoothParam,maxNiter)
    
    if nargin < 6
        maxNiter = 50;
    end
    
    if size(mbdata,2) ~= sum(sum(fExists))
        error('error: WTF Seriously.')
    end
    if sum(sum((fExists+emptySpec)==1)) ~= sum(sum(fExists)) + sum(sum(emptySpec))
        error('error: Seriously. Fucking Seriously.')
    end
    
    fExists_new = (fExists+emptySpec) == 1;
    % Create XYL
    XYL = double(fExists_new);
    XYL(fExists_new) = 1:sum(sum(fExists_new));
    % Remember which ones are the old spectra
    spec_old = XYL(fExists);
    fExists = fExists_new;
    clear fExists_new
    % Add empty spectra to mbdata
    temp = mbdata;
    mbdata = false(size(mbdata,1),sum(sum(fExists)));
    mbdata(:,spec_old) = temp;
    clear temp spec_old
    
	% Create LXY
	LXY = zeros(sum(sum(fExists)),3);
    for spec_idx = 1:size(LXY,1)
        [i,j] = ind2sub(size(XYL),find(XYL == spec_idx));
        LXY(spec_idx,1) = i;
        LXY(spec_idx,2) = j;
    end

	% Extend to include empty border
	XYLext = zeros(size(XYL,1)+2,size(XYL,2)+2);
    XYLext(2:1+size(XYL,1),2:1+size(XYL,2)) = XYL;
	
	% Begin Smoothing
	continueSmoothing = true(1,1);
    nSmooths = 0;
    while continueSmoothing
		nSmooths = nSmooths + 1;
		temp = mbdata;

        mNsame = zeros(size(mbdata,1),size(mbdata,2));
		mNadja = zeros(1,size(mbdata,2));
        for spec = 1:size(mbdata,2)
			x_cur = LXY(spec,1);
			x_cur = x_cur + 1;
			y_cur = LXY(spec,2);
			y_cur = y_cur + 1;
			XYLzoom = XYLext((x_cur-1):(x_cur+1),(y_cur-1):(y_cur+1));
			specOfInterest = XYLzoom(XYLzoom > 0);
			mNadja(spec) = length(specOfInterest);         
        %             disp(['      Spectrum ' num2str(spec) ' - ' num2str(mNadja(spec)) ' adjacent spectra'])            
			curVal = mbdata(:,spec);
			mNsame(:,spec) = (1-curVal).*repmat(mNadja(spec),length(curVal),1) + (2*curVal - 1).*sum(mbdata(:,specOfInterest),2);
        end
        uniMeasure9vals = zeros(size(mbdata));
        uniMeasure9vals(:,mNadja > 1) = (mNsame(:,mNadja > 1) - 1)./repmat(mNadja(mNadja > 1)-1,size(mbdata,1),1);
					
		mbdata(uniMeasure9vals <= smoothParam) = 1 - mbdata(uniMeasure9vals <= smoothParam);
%         disp([num2str(nSmooths) '. Number of changes: ' num2str(sum(sum(uniMeasure9vals <= smoothParam))) ' (' num2str(sum(sum(uniMeasure9vals <= smoothParam))*100/size(mbdata,1)/size(mbdata,2)) '%)'])
% 		disp(['     n: ' num2str(size(mbdata,2) - sum(sum(mbdata) == 0)) '   d: ' num2str(size(mbdata,1) - sum(sum(mbdata,2) == 0))])

		if sum(sum(mbdata ~= temp)) == 0
%                     disp(' converged.')
            continueSmoothing = false(1,1);
            converged = true(1,1);
		elseif nSmooths == maxNiter
% 			disp(' ---------------------------------- ')
			disp(['warning: smooth did not converge in ' num2str(nSmooths) ' iterations.'])
% 			disp(' ---------------------------------- ')
			continueSmoothing = false(1,1);
            converged = false(1,1);
		end

    end
    
    nIter = nSmooths;
    disp(['  ' num2str(smoothParam) '-smooth converged after ' num2str(nSmooths) ' iterations.'])
	
	emptyObs = sum(mbdata) == 0;
	emptySpec = fExists;
	emptySpec(fExists) = emptyObs;
	fExists(fExists) = ~emptyObs;
	emptyDims = sum(mbdata,2) == 0;
	vbincentrs = vbincentrs(~emptyDims);
	mbdata = mbdata(~emptyDims,~emptyObs);
	
end














