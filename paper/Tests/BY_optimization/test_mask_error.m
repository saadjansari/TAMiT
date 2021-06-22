function simImg = test_mask_error( imageIn, spindleObj)


    simImg=1;
    
    curve = spindleObj.featureList{3}.featureList{2};
    
    % Increase curve length until it gets inside mask
    cc = round( curve.GetCoords() );
    idx = sub2ind( size( spindleObj.image), cc(2,:), cc(1,:), cc(3,:) );
    imgDummy = zeros( size(spindleObj.image) );
    imgDummy( idx) = 1;
    m_invert = ~spindleObj.mask;
    
    outsideXY = sum( imgDummy(:) .* m_invert(:) ) > 0;
    while ~outsideXY
        
        curve.L = 1.05*curve.L;
        cc = round( curve.GetCoords() );
        idx = sub2ind( size( spindleObj.image), cc(2,:), cc(1,:), cc(3,:) );
        imgDummy = zeros( size(spindleObj.image) );
        imgDummy( idx) = 1;
        outsideXY = sum( imgDummy(:) .* m_invert(:)) > 0;
    end
    
    figure; imagesc( max(m_invert.*imgDummy, [],3))
    
    % Now implement error routine for this curve
    imgOutside = logical( m_invert.*imgDummy);
    imgOutsideDilated = max( imdilate( imgOutside, strel('disk',4) ), [], 3);
    
    % Get indices that lie inside mask
    idx = find(imgOutsideDilated);
    [yy,xx] = ind2sub( size(imgOutsideDilated), idx);
    
    % simulated curve image
    simImg = curve.simulateFeature( size(imageIn) ) .* ;
    
    
    
    
    err=1;
    % find boundary of 2Dmask
    
    figure; imagesc( max(m_invert,[],3)); hold on;
    plot( B{1}(:,2), B{1}(:,1), 'color', 'w', 'linewidth',3)
    

end