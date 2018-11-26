function stmt = convertVec2Struct_MT( vec, pStruct)
% converts a vector output, possibly the final fitted vector through lsqnonlin, into a structure where the vector is parsed accurately

try
    nummt = pStruct.numberOfMicrotubules;
    dimmt = pStruct.dim;
    fixmtoc = pStruct.StartingPointFixed;
    pOrder = pStruct.PolyOrder;
    if dimmt==3, pOrderZ = pStruct.PolyOrderZ; end
catch, error( 'input parameter structure is insufficient for calculation'), end

if dimmt ~= 2 && dimmt ~= 3, error('mt must exist in a 2D or 3D image'), end

% figure out the number of parameters for a single microtubule
numParMT = (length(vec)-1 ) / nummt; 
if mod( numParMT, 1) ~= 0, error( 'size of input vector is incorrect'), end

if dimmt==2
    numParMT_also = 1 + dimmt + dimmt*(pOrder+1); if fixmtoc, numParMT_also = numParMT_also-dimmt; end
elseif dimmt==3
    numParMT_also = 1 + 3 + 2*(pOrder+1) + 1*(pOrderZ+1); if fixmtoc, numParMT_also = numParMT_also-dimmt; end
end


if numParMT ~= numParMT_also, error('Issue with the parameters supplied in the struct. MT param vec length is inconsistent with them.'), end

% figure out number of parameters for coef per dimension of a microtubule
if fixmtoc
    ncoefXY = pOrder+1-1; 
    if dimmt==3, ncoefZ = pOrderZ+1-1; end 
else 
    ncoefXY = pOrder+1; 
    if dimmt==3, ncoefZ = pOrderZ+1; end 
end

stmt = pStruct;
stmt.background = vec(1);
for jmt = 1 : nummt
    cid = 2 + numParMT*(jmt-1); initid = cid;
    stmt.amplitude(jmt) = vec( cid); cid=cid+1; 
    stmt.std{jmt} = vec( cid:cid+dimmt-1); cid = cid+dimmt;
    stmt.coefX{jmt} = vec( cid:cid+ncoefXY-1); cid = cid+ncoefXY;
    stmt.coefY{jmt} = vec( cid:cid+ncoefXY-1); cid = cid+ncoefXY;
    if dimmt==3, stmt.coefZ{jmt} = vec( cid:cid+ncoefZ-1); cid = cid+ncoefZ; end
    if fixmtoc, 
        stmt.coefX{jmt} = [stmt.coefX{jmt}, pStruct.XCoefEnd];
        stmt.coefY{jmt} = [stmt.coefY{jmt}, pStruct.YCoefEnd];
        if dimmt==3, stmt.coefZ{jmt} = [stmt.coefZ{jmt}, pStruct.ZCoefEnd]; end
    end
    if numParMT ~= (cid-initid), error('Issue with propagation of current index in for loop of value assignment to structure.'), end
end

stmt.numberOfMicrotubules = nummt;
stmt.dimmt = dimmt;
stmt.polyOrder = pOrder;
if dimmt==3, stmt.polyOrderZ = pOrderZ;
stmt.startPointFixed = fixmtoc;

end
