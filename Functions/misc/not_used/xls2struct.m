function st  = xls2struct( xlsfile)
%  xls2struct: assigns each column in an xls file to the field given by column headers in the xls file.

[ data, headers] = xlsread( xlsfile);

headers = headers(1,:);
numHeaders = length( headers);

% parse the labels into structure var names and values into their place
for jHead = 1 : numheaders
    st.headers{ jHead } = data( :, jHead);
end

end
