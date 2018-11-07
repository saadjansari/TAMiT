function stmt = convertVec2Struct_MT( vec, pStruct)
% converts a vector output, possibly the final fitted vector through lsqnonlin, into a structure where the vector is parsed accurately

try
    nummt = pStruct.numberOfMicrotubules;
    dimmt = pStruct.dim;
    fixmtoc = pStruct.StartingPointFixed;
    pOrder = pStruct.PolyOrder;
catch, error( 'input parameter structure is insufficient for calculation'), end

if dimmt ~= 2 && dimmt ~= 3, error('mt must exist in a 2D or 3D image'), end

% figure out the number of parameters for a single microtubule
numParMT = (length(vec)-1 ) / nummt; 
if mod( numParMT, 1) ~= 0, error( 'size of input vector is incorrect'), end
numParMT_also = 1 + dim + dim*(pOrder+1); if fixmtoc, numParMT_also = numParMT_also-1; end
if numParMT ~= numParMT_also, error('Issue with the parameters supplied in the struct. MT param vec length is inconsistent with them.'), end

% figure out number of parameters for coef per dimension of a microtubule
if fixmtoc, ncoefmt = pOrder+1-1; else ncoefmt = pOrder+1; end

stmt = pStruct;
stmt.background = vec(1);
for jmt = 1 : nummt
    cid = 1 + numParMT*(jmt-1); initid = cid;
    stmt.amplitude(jmt) = vec( cid); cid=cid+1; 
    stmt.std{jmt} = vec( cid:cid+dim-1); cid = cid+dim;
    stdmt.coefX{jmt} = vec( cid:cid+ncoefmt-1); cid = cid+ncoefmt;
    stdmt.coefY{jmt} = vec( cid:cid+ncoefmt-1); cid = cid+ncoefmt;
    stdmt.coefZ{jmt} = vec( cid:cid+ncoefmt-1); cid = cid+ncoefmt;
    if fixmtoc, 
        stdmt.coefX{jmt} = [stdmt.coefX{jmt}, pStruct.XCoefEnd];
        stdmt.coefY{jmt} = [stdmt.coefY{jmt}, pStruct.YCoefEnd];
        stdmt.coefZ{jmt} = [stdmt.coefZ{jmt}, pStruct.ZCoefEnd];
    end





    % Add condition for 3rd dimension






    if numParMT ~= (cid-initid-1), error('Issue with propagation of current index in for loop of value assignment to structure.'), end
end


