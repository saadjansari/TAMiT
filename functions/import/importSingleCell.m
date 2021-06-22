function cellData = importSingleCell( filepath)
  
%{
Ensure that the file is of a valid type
%}
    validFormats = {'.mat','.tiff'};
    [ rpath, fname, ext] = fileparts( filepath);
    if ~any( strcmp( ext, validFormats) )
        error('importSingleCell : file extension is invalid'); end

    % Import depending on file extension
    switch ext
        case '.mat'
            %{
            Expected data structure in filepath
            CellData : main variable
            CellData.raw (3D raw data)
            CellData.cell3D (3D masked data) (optional)
            CellData.planeTimes (average time values for all XY slices) (optional)
            CellData.cellNumber (cell number in segmentations routine) (optional)
            CellData.locations (image with cell numbers and locations shown) (optional)
            CellData.cellCentroids (array with XY positions of cells masked in original segmentation routine) (optional)
            CellData.metaData (metaData)
            %}
            load(filepath);
        case '.tiff'
            error('importSingleCell : .tiff not coded')
    end


end
