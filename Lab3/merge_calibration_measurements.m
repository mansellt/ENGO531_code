function merge_calibration_measurements()
% Combine independently measured project CSVs into one calibration.PHO.
% Each teammate works in a SEPARATE project folder. Do not edit a shared CSV
% simultaneously. All photo files must be identical and named Image1.jpg etc.

[files, folder] = uigetfile('*.csv', ...
    'Select calibration_points.csv files from each collaborator', ...
    'MultiSelect', 'on');
if isequal(files,0), return; end
if ischar(files), files = {files}; end
% uigetfile MultiSelect normally selects files in ONE directory. For files in
% different directories, collect additional file paths interactively.
paths = cellfun(@(f) fullfile(folder,f),files,'UniformOutput',false);
while true
    answer = questdlg('Add a calibration_points.csv from another folder?', ...
        'More collaborators','Add file','Done','Done');
    if ~strcmp(answer,'Add file'),break;end
    [f,p] = uigetfile('*.csv','Select another collaborator CSV');
    if isequal(f,0),break;end
    paths{end+1}=fullfile(p,f); %#ok<AGROW>
end
T = table(strings(0,1),zeros(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'point_id','image_id','filename','x_pix','y_pix'});
for k = 1:numel(paths)
    A = readtable(paths{k},'Delimiter',',','FileType','text', ...
        'TextType','string','VariableNamingRule','preserve');
    A.Properties.VariableNames=strtrim(A.Properties.VariableNames);
    required={'point_id','image_id','filename','x_pix','y_pix'};
    assert(all(ismember(required,A.Properties.VariableNames)), ...
        'Bad columns in %s',paths{k});
    A=A(:,required);
    A.point_id=string(A.point_id);
    A.image_id=double(A.image_id);
    A.filename=string(A.filename);
    A.x_pix=double(A.x_pix); A.y_pix=double(A.y_pix);
    for j=1:height(A)
        [~,stem,ext]=fileparts(char(A.filename(j)));
        token=regexp(stem,'^Image(\d+)$','tokens','once','ignorecase');
        assert(~isempty(token) && str2double(token{1})==A.image_id(j), ...
            'Image ID mismatch in %s row %d',paths{k},j);
        A.filename(j)=string([stem ext]);
    end
    T=[T;A]; %#ok<AGROW>
end
assert(~isempty(T),'No measurements found.');
keys = string(T.image_id) + "|" + T.point_id;
[u,~,group]=unique(keys);
count=accumarray(group,1);
if any(count>1)
    bad=u(count>1);
    error(['Duplicate measurements for image|point combinations: %s. ' ...
        'Assign different images to different people or resolve conflicts before merging.'], ...
        strjoin(cellstr(bad),', '));
end
T=sortrows(T,{'image_id','point_id'});
outputFolder=uigetdir(pwd,'Choose output folder for MERGED measurements');
if isequal(outputFolder,0),return;end
outCSV=fullfile(outputFolder,'calibration_points.csv');
outPHO=fullfile(outputFolder,'calibration.PHO');
outLookup=fullfile(outputFolder,'image_id_lookup.csv');
if isfile(outCSV)||isfile(outPHO)||isfile(outLookup)
    answer=questdlg('Output files already exist. Overwrite them?', ...
        'Confirm overwrite','Overwrite','Cancel','Cancel');
    if ~strcmp(answer,'Overwrite'),return;end
end
L=unique(T(:,{'image_id','filename'}),'rows');
assert(numel(unique(L.image_id))==height(L), ...
    'Multiple filenames have the same image ID. Resolve first.');
writetable(T,outCSV);
writetable(L,outLookup);
[fid,msg]=fopen(outPHO,'w');
assert(fid>=0,'Cannot open PHO: %s',msg);
cleanup=onCleanup(@() fclose(fid));
for j=1:height(T)
    fprintf(fid,'%s %d %.2f %.2f\n',char(T.point_id(j)), ...
        T.image_id(j),T.x_pix(j),T.y_pix(j));
end
fprintf('Merged %d observations from %d CSV files for %d images.\n', ...
    height(T),numel(paths),height(L));
end
