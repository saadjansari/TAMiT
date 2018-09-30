function coordsMT  = estimateMicrotubulePoints( startPt, startOrient, imageMT, stepSize, visibility, fieldOfVision)
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

            % clear out any previously stored estimated points
            
            success = 1; iter = 1;
            max_mt_length = 80; % max MT length approx
            max_iter = round(max_mt_length/stepSize);
            orientationOld = startOrient;
            coordCurr = startPt;
            coordsMT = coordCurr; 
            orientBank = [orientationOld];
            % iteratively propagate along microtubule
            while success && iter< max_iter
                
                if iter <= 2, fov = 1.5*fieldOfVision; else, fov = fieldOfVision; end

                [coordNext, success, orientationNext] = estimateNextPoint( coordCurr, orientationOld, stepSize, visibility, fov, imageMT);
                if success
%                     plotEstimatedPointNext( obj, orientationNext, stepSize, visibility, fieldOfVision, iter)
%                     drawnow; pause(0.2)
                    orientBank = [ orientBank, orientationNext];
                    try, orientationOld = orientBank(end-2:end); end, try, orientationOld = orientBank(end-1:end); end
                    coordCurr = coordNext;
                    coordsMT = [ coordsMT , coordNext];

                end
                
                iter = iter+1;
                
            end
            
            % try continuing with increased visibility (when there is an lateral overlap of two MTS, bigger visibility can help decide the best path 
            success = 1; count=iter;
            
            while success && iter < count+2 
                
                if iter <= 2, fov = 1.5*fieldOfVision; else, fov = fieldOfVision; end

                [coordNext, success, orientationNext] = estimateNextPoint( coordCurr, orientationOld, stepSize, visibility*1.5, fov, imageMT);
            
                if success
%                     plotEstimatedPointNext( obj, orientationNext, stepSize, visibility, fieldOfVision, iter)
%                     drawnow; pause(0.2)
                    orientBank = [ orientBank, orientationNext];
                    try, orientationOld = orientBank(end-2:end); end, try, orientationOld = orientBank(end-1:end); end
                    coordCurr = coordNext;
                    coordsMT = [ coordsMT, coordNext];
                end
                iter = iter+1;
                
            end
            
            
            % try continuing with half the stepSize, fine ending (ensuring the maximum possible length is reached)
            success = 1; count=iter;
            
            while success && iter < count+2 
                
                if iter <= 2, fov = 1.5*fieldOfVision; else, fov = fieldOfVision; end

                [coordNext, success, orientationNext] = estimateNextPoint( coordCurr, orientationOld, stepSize, visibility/2, fov, imageMT);
            
                if success
%                     plotEstimatedPointNext( obj, orientationNext, stepSize, visibility, fieldOfVision, iter)
%                     drawnow; pause(0.2)
                    orientBank = [ orientBank, orientationNext];
                    try, orientationOld = orientBank(end-2:end); end, try, orientationOld = orientBank(end-1:end); end
                    coordCurr = coordNext;
                    coordsMT = [ coordsMT, coordNext];
                end
                iter = iter+1;
                
            end
end
