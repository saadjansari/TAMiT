function [frames_padded, feats_padded] = padding_on( frames, feats, pad_xy, pad_z)

    feats_padded = feats;
    frames_padded = frames;
    for jfeat = 1 : length(feats_padded)
        
        feat = feats_padded{jfeat}.copyDeep();
        frame = frames{jfeat};
        % Pad image
        img0 = frame;
        framed = zeros( 2*pad_xy+size(img0, 1), 2*pad_xy+size(img0,2), 2*pad_z+size(img0,3), class(img0) );
        framed( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z) = img0;

        % Pad features
        % Main organizer
        feat.image = framed;
        mask0 = feat.mask;
        mask1 = zeros( size(framed,1), size(framed,2), size(framed,3), class(mask0) );
        mask1( 1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1+pad_z:end-pad_z) = mask0;
        mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, 1) = mask0(:,:,1);
        mask1(1+pad_xy:end-pad_xy, 1+pad_xy:end-pad_xy, end) = mask0(:,:,end);
        feat.mask = mask1;

        fspindle = feat.featureList{1};
        fspindle.startPosition = fspindle.startPosition + [pad_xy,pad_xy,pad_z];
        fspindle.endPosition = fspindle.endPosition + [pad_xy,pad_xy,pad_z];
        fspindle.bounds.ub.startPosition = fspindle.bounds.ub.startPosition + 2*[pad_xy,pad_xy,pad_z];
        fspindle.bounds.ub.endPosition = fspindle.bounds.ub.endPosition + 2*[pad_xy,pad_xy,pad_z];

        for ja = 2 : length( feat.featureList)
            faster =  feat.featureList{ja};
            for jf = 1 : faster.numFeatures

                if strcmp( faster.featureList{jf}.type, 'Spot')
                    faster.featureList{jf}.position = faster.featureList{jf}.position + [pad_xy,pad_xy,pad_z];
                    faster.featureList{jf}.bounds.ub.position = faster.featureList{jf}.bounds.ub.position + 2*[pad_xy,pad_xy,pad_z];
                elseif strcmp( faster.featureList{jf}.type, 'CurvedMT')
                    faster.featureList{jf}.origin = faster.featureList{jf}.origin' + [pad_xy,pad_xy,pad_z];
                    faster.featureList{jf}.bounds.ub.origin = faster.featureList{jf}.bounds.ub.origin + 2*[pad_xy,pad_xy,pad_z];
                else
                    error('unknown feature type')
                end
                faster.featureList{jf}.fillParams( size(framed) );
            end
        end
        frames_padded{ jfeat} = framed;
        feats_padded{jfeat} = feat;
    end
end
