function [vec, vecUb, vecLb] = convertStruct2Vec_MT( stmt, stmt_ub, stmt_lb, fixedMTOC, fixCoef)
%  converts a structure with information about the features and its environment to a vector that can be fed into MATLAB's lsqnonlin toolbox

%  the vector looks like this:
%  vec = [ background, mt1, mt2, ... ]
%  mtX = [ amplitude, stdX, stdY, stdZ, CoefX, CoefY, CoefZ ];

if nargin < 5
    fixXY = 0;
elseif strcmp( fixCoef, 'fixXY')
    fixXY = 1;
else
    fixXY = 0;
end

vec = stmt.background;
vecUb = stmt_ub.background;
vecLb = stmt_lb.background;

if fixedMTOC, cRge = 1:length(stmt.coefX{1})-1; else cRge = 1:length(stmt.coefX{1}); end
if stmt.dim==3 && fixedMTOC, cRgeZ = 1:length(stmt.coefZ{1})-1; elseif stmt.dim==3 && ~fixedMTOC, cRgeZ = 1:length(stmt.coefZ{1}); end

stmt.dim
switch stmt.dim
    case 2
        for jmt = 1 : stmt.numberOfMicrotubules
            vec = [vec, stmt.amplitude(jmt), stmt.std{jmt}, stmt.coefX{jmt}(cRge), stmt.coefY{jmt}(cRge) ];
            vecUb = [vecUb, stmt_ub.amplitude(jmt), stmt_ub.std{jmt}, stmt_ub.coefX{jmt}(cRge), stmt_ub.coefY{jmt}(cRge)];
            vecLb = [vecLb, stmt_lb.amplitude(jmt), stmt_lb.std{jmt}, stmt_lb.coefX{jmt}(cRge), stmt_lb.coefY{jmt}(cRge)];
        end 
    case 3
        for jmt = 1 : stmt.numberOfMicrotubules
            if fixXY
                vec = [vec, stmt.amplitude(jmt), stmt.std{jmt}, stmt.coefX{jmt}(cRge), stmt.coefY{jmt}(cRge), stmt.coefZ{jmt}(cRgeZ) ];
                vecUb = [vecUb, stmt_ub.amplitude(jmt), stmt_ub.std{jmt}, stmt.coefX{jmt}(cRge), stmt.coefY{jmt}(cRge), stmt_ub.coefZ{jmt}(cRgeZ) ];
                vecLb = [vecLb, stmt_lb.amplitude(jmt), stmt_lb.std{jmt}, stmt.coefX{jmt}(cRge), stmt.coefY{jmt}(cRge), stmt_lb.coefZ{jmt}(cRgeZ) ];
            else
                vec = [vec, stmt.amplitude(jmt), stmt.std{jmt}, stmt.coefX{jmt}(cRge), stmt.coefY{jmt}(cRge), stmt.coefZ{jmt}(cRgeZ) ];
                vecUb = [vecUb, stmt_ub.amplitude(jmt), stmt_ub.std{jmt}, stmt_ub.coefX{jmt}(cRge), stmt_ub.coefY{jmt}(cRge), stmt_ub.coefZ{jmt}(cRgeZ) ];
                vecLb = [vecLb, stmt_lb.amplitude(jmt), stmt_lb.std{jmt}, stmt_lb.coefX{jmt}(cRge), stmt_lb.coefY{jmt}(cRge), stmt_lb.coefZ{jmt}(cRgeZ) ];
            end
        end 
    otherwise
        error('Microtubules must exist in either 2D or 3D images')
end

end
