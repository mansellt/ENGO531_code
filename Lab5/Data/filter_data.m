%
data = readtable("ENGO531_Lab5_Group3.txt");

data(:,4:10) = [];

writetable(data, "filtered_data");