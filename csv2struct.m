function st  = csv2struct( csvfile)

dat = importdata( csvfile);

labels = dat.colheaders;
values = dat.data;
numLabels = length( labels);
iterVals = size(values, 1);

% parse the labels into structure var names and values into their place
for jLab = 1 : numLabels
    st.labels{jLab} = values( :, jLab);
end

end
