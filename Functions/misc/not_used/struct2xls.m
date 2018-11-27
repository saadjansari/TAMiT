function st  = data2xls( xlsName, savePath, headers, data)
% data2xls: Saves headers into the first row of the xls file, and save data in rows 2 onward. Each coilumn of data corresponds to an entry in headers

xlsPath = [savePath, filesep, xlsName];
xlswrite( xlsPath , [headers; data] );

end
