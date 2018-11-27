function D=createDistanceMatrixPolar(M,N)
% createDistanceMatrixPolar calculates the distance matrix for two sets of
% points in polar or spherical polar coordinates
%
% SYNOPSIS   D=createDistanceMatrixPolar(M,N)
%
% INPUT      M and N are the matrices containing the set of polar point coordinates.
%            M and N can represent point positions in 1, 2 and 3D, as follows.
%            
%            In 1D: M=[ R1        and   N=[ R1
%                       R2                  R2
%                       ...                 ... 
%                       Rm ]                Rn ]
%
%            In 2D:
%                   M=[ R1 th1     and   N=[ R1 th1
%                       R2 th2               R2 th2
%                        ...                 ...
%                       Rm thm ]             Rn thn ]
%
%            In 3D:
%                   M=[ R1 th1 phi1  and   N=[ R1 th1 phi1
%                       R2 th2 phi2            R2 th2 phi2
%                         ...                  ...
%                       Rm thm phim ]          Rn thn phin ]
%
%
% OUTPUT   D : distance matrix D=(dij), i=1..m, j=1..n
% 
% REMARK   For 1D, both positive and negative distances are returned.
%
% Meredith Betterton 12/22/15

%get dimension of problem
dimension=size(M,2);

%if two vectors have different dimensions, return error
if size(N,1) ~= dim
    error('createDistanceMatrixPolar: different dimensions of input vectors');
end

%get lengths of each vector
m=size(M,1);
n=size(N,1);

switch dimension
    case 1 %1D problem, take differences of lengths
        Mext=repmat(M,1,n);
        Next=repmat(N',m,1);
        D=Mext-Next;
    case 2
        %convert R, theta to x, y: x=R cos theta, y=R sin theta
        mXcoord=repmat(M(1,:).*cos(M(2,:)),1,n);
        mYcoord=repmat(M(1,:).*sin(M(2,:)),1,n);
        nXcoord=repmat((N(1,:).*cos(N(2,:)))',m,1);
        nYcoord=repmat((N(1,:).*sin(N(2,:)))',m,1);
        %calculate distance
        D=sqrt( (mXcoord-nXcoord).^2 + (mYcoord-nYcoord).^2 );
    case 3
        %convert R, theta, phi to x, y, z: x=R sin theta cos phi, y=R sin
        %theta sin phi, z=R cos theta
        mXcoord=repmat(M(1,:).*sin(M(2,:)).*cos(M(3,:)),1,n);
        mYcoord=repmat(M(1,:).*sin(M(2,:)).*sin(M(3,:)),1,n);
        mZcoord=repmat(M(1,:).*cos(M(2,:)),1,n);
        nXcoord=repmat((N(1,:).*sin(N(2,:)).*cos(N(3,:)))',m,1);
        nYcoord=repmat((N(1,:).*sin(N(2,:)).*sin(N(3,:)))',m,1);
        nZcoord=repmat((N(1,:).*cos(N(2,:)))',m,1);
        %calculate distance
        D=sqrt( (mXcoord-nXcoord).^2 + (mYcoord-nYcoord).^2 + (mZcoord-nZcoord).^2 );
    otherwise
        error('createDistanceMatrixPolar: wrong dimensions of input vectors');
end
    