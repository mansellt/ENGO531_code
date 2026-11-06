%
data = readtable("ENGO531_Lab5_Group3_New.txt");

data(:,4:10) = [];

writetable(data, "filtered_data");