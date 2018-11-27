function varargout = joinColumns(delimeter,varargin)
%joinColumns vertically concatenate columns separated by a delimeter
%
% [colVecA, colVecB, colVecC, ...] = joinColumns(delimeter, A, B, C, ...)
%
% INPUT
%
% delimeter is a scalar value that will be placed between the columns of
%    the matrices. An empty matrix indicates no delimeter.
%    (optional. A non-scalar value is considered a matrix and the delimeter
%    will be empty)
%
% A, B, C ... are matrices that can hold delimeter as a value. They need
% not be the same size or type
%
% OUTPUT
%
% columnVector - a column vector is returned for each matrix generated by
% concatenating the matrices vertically separated by the scalar delimeter
%
% joinColumns(matrix) and joinColumns([],matrix) are equivalent to matrix(:)
    
    % INPUT VALIDATION: If delimeter is not a scalar value, then assume no
    % delimeter is given and that the matrix should be turned into a column
    % vector
    if(~isscalar(delimeter))
        % add delimeter to the cell array of matrices to linearize
        % if delimeter is empty, this should have no effect
        varargin = [ delimeter varargin];
        delimeter = [];
    end
    
    % For each matrix given linearize and output. This allows for
    % vectorized syntax and is compatible with paralleization with
    % distributed objects
    
    % We only need to check if the delimeter is empty once
    if(~isempty(delimeter))
        varargin = cellfun(@appendDelimeter,varargin,'UniformOutput',false);
        varargout = cellfun(@joinColumnsSingle,varargin,'UniformOutput',false);
        varargout = cellfun(@(x) x(1:end-1),varargout,'UniformOutput',false);
    else
        varargout = cellfun(@joinColumnsSingle,varargin,'UniformOutput',false);
    end
    
    function matrix = appendDelimeter(matrix)
        % append the delimeter as the last row
        matrix(end+1,:) = delimeter;
    end
    function columnVector = joinColumnsSingle(matrix)
        % linearize
        columnVector = matrix(:);
    end
end