function coordsMT  = estimateMicrotubulePoints( startPt, startOrient, imageMT, stepSize, visibility, fieldOfVision, bkg_thresholding, plotflags)
            % EstimateMicrotubulePoints: estimates points along
            % microtubule.
            % We take the following approach:
            % 1) Nucleate at the MTOC. We find the angle that gives the maximum in a
            %    radially integrated signal, up to some max radius( which we will call
            %    the visibility)
            % 3) We'll propagate from the MTOC to a new point along the optimum angle,
            %    which a defined distance away from the MTOC (this we will call the
            %    step size). Throughout this process, we might look for optimum angles
            %    within a certain range related to the orientation of the Tube (this we
            %    will refer to as the field of vision).
            % 4) We will continue propagating until we can't find an optimum angle
            %    anymore (all our surroudings appear the same. This will be based on
            %    some optimality condition. At this point we will change our step size
            %    to half its original value and see if we can propogate a smaller
            %    distance forward.
            % 5) This will conclude our initial estimation phase.

            if nargin < 7, bkg_thresholding = 1; end
            if nargin < 8, plotflags.success = 0; plotflags.fail=0; end
            
            % clear out any previously stored estimated points
             
            success = 1; iter = 1;
            max_mt_length = 100; % max MT length approx
            max_iter = round(max_mt_length/stepSize);
            orientationOld = startOrient;
            coordCurr = startPt;
            coordsMT = coordCurr; 
            orientBank = [orientationOld];
            totalImprovedCones = 2;
            currImproveCones = 0;
            % iteratively propagate along microtubule
            while success && iter< max_iter && currImproveCones <= totalImprovedCones
                
                if iter <= 2, fov = 1.5*fieldOfVision; else, fov = fieldOfVision; end

                [coordNext, success, orientationNext] = estimateNextPoint( coordCurr, orientationOld, stepSize, visibility, fov, imageMT, bkg_thresholding, plotflags);
                if ~success % try again for a limited number of times with an improved search cone
                    [coordNext, success, orientationNext] = estimateNextPoint( coordCurr, orientationOld, stepSize, visibility/2, 1.5*fov, imageMT, bkg_thresholding, plotflags);
                    currImproveCones = currImproveCones+1;
                end

                if success
                    orientBank = [ orientBank, orientationNext];
                    try, orientationOld = orientBank(end-2:end); end, try, orientationOld = orientBank(end-1:end); end
                    coordCurr = coordNext;
                    coordsMT = [ coordsMT , coordNext];

                end

                iter = iter+1;
                
            end
            
end
