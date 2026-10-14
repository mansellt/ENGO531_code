function measure_image_points_resumable_v2()
% ENGO 531 Lab 3: resumable manual image-point measurements.
% Run: measure_image_points_resumable_v2
%
% One image at a time; measurements saved after EVERY add/edit/delete.
% You can quit MATLAB and resume later by choosing the SAME project folder.
%
% Output files in project folder:
%  calibration.PHO         point_id image_id x y (no header; two decimals)
%  calibration_points.csv  editable record of measurements
%  image_id_lookup.csv     image IDs and paths (keep this with your data)
%  image_XX_annotated.png  inspection image (refreshed when closing image)
%
% Pixel system: upper-left origin, x right, y downward. Bundle INT y_dir=-1.
% Use unrotated, uncropped, full-resolution photographs consistently.
%
% Requires MATLAB's image display and drawpoint tools.

projectFolder = uigetdir(pwd,'Select a folder to KEEP your calibration measurements');
if isequal(projectFolder,0), return; end
pointsPath = fullfile(projectFolder,'calibration_points.csv');
lookupPath = fullfile(projectFolder,'image_id_lookup.csv');

if isfile(pointsPath)
    T = readtable(pointsPath,'Delimiter',',','FileType','text', ...
        'TextType','string','VariableNamingRule','preserve');
    T.Properties.VariableNames = strtrim(T.Properties.VariableNames);
    requiredT = {'point_id','image_id','filename','x_pix','y_pix'};
    assert(all(ismember(requiredT,T.Properties.VariableNames)), ...
        'Unexpected calibration_points.csv columns: %s',strjoin(T.Properties.VariableNames,', '));
    T.point_id = string(T.point_id);
    T.filename = string(T.filename);
    T.image_id = double(T.image_id);
    T.x_pix = double(T.x_pix);
    T.y_pix = double(T.y_pix);
else
    T = table(strings(0,1),zeros(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'point_id','image_id','filename','x_pix','y_pix'});
end
if isfile(lookupPath)
    L = readtable(lookupPath,'Delimiter',',','FileType','text', ...
        'TextType','string','VariableNamingRule','preserve');
    L.Properties.VariableNames = strtrim(L.Properties.VariableNames);
    requiredL = {'image_id','filename'};
    assert(all(ismember(requiredL,L.Properties.VariableNames)), ...
        'Unexpected image_id_lookup.csv columns: %s',strjoin(L.Properties.VariableNames,', '));
    L.filename = string(L.filename);
    L.image_id = double(L.image_id);
else
    L = table(zeros(0,1),strings(0,1), ...
        'VariableNames',{'image_id','filename'});
end

% Do not overwrite any existing files just from opening a project.
% Existing data will be saved after an add, edit, delete, or new image registration.
fprintf('Loaded %d observations and %d image registrations from %s\n', ...
    height(T),height(L),projectFolder);
% Select an image folder, rather than using macOS filtered file chooser.
% This avoids the issue where PNG/JPEG photographs appear greyed out.
imageFolder = uigetdir(projectFolder, ...
    'Select the FOLDER containing your PNG or JPEG photographs');
if isequal(imageFolder,0), return; end

while true
    % List valid photographs in MATLAB itself. No file-type dialog filters.
    contents = dir(imageFolder);
    contents = contents(~[contents.isdir]);
    names = {contents.name};
    valid = ~cellfun('isempty',regexpi(names,'\.(png|jpg|jpeg|tif|tiff|bmp)$'));
    names = names(valid);
    if isempty(names)
        answer=questdlg('No PNG/JPEG/TIFF/BMP files were found in this folder.', ...
            'No images','Choose another folder','Quit','Choose another folder');
        if ~strcmp(answer,'Choose another folder'),break;end
        imageFolder = uigetdir(imageFolder,'Select the image folder');
        if isequal(imageFolder,0),break;end
        continue;
    end
    names=sort(names);
    options=[names, {'[CHOOSE A DIFFERENT FOLDER]','[FINISH]'}];
    [selection,ok]=listdlg('Name','Choose a photograph', ...
        'PromptString','Select ONE image (or change folder):', ...
        'ListString',options,'SelectionMode','single', ...
        'ListSize',[480 360]);
    if ~ok || selection==numel(options),break;end
    if selection==numel(options)-1
        newFolder=uigetdir(imageFolder,'Select the image folder');
        if ~isequal(newFolder,0),imageFolder=newFolder;end
        continue;
    end
    filename=names{selection};
    imagePath=string(fullfile(imageFolder,filename));
    % Read first; register the photograph only if decoding succeeded.
    try
        im = imread(char(imagePath));
    catch ME
        warning('Cannot open %s: %s',imagePath,ME.message);
        warndlg(sprintf('Could not read this image:\n%s',ME.message), ...
                'Image read error');
        continue;
    end
    ix = find(L.filename==imagePath,1);
    if isempty(ix)
        if isempty(L), idImage=1; else, idImage=max(L.image_id)+1; end
        L(end+1,:) = {idImage,imagePath};
        saveData();
    else
        idImage=L.image_id(ix);
    end

    info = imfinfo(char(imagePath));
    if isfield(info,'Orientation') && isnumeric(info.Orientation) && ...
            ~isempty(info.Orientation) && info.Orientation~=1
        warning('EXIF orientation is %d. Verify orientation is correct.',info.Orientation);
    end
    fig = figure('Name',sprintf('Image %d: %s',idImage,filename), ...
        'NumberTitle','off','Color','w');
    ax = axes('Parent',fig);
    imshow(im,'Parent',ax,'InitialMagnification','fit');
    hold(ax,'on');
    refreshLabels();

    while isgraphics(fig)
        choice = menu(sprintf('Image %d: %s | %d measurements', ...
            idImage,filename,sum(T.image_id==idImage)), ...
            'Add a point','Edit an existing point','Delete a point', ...
            'Save and choose another image','Save and quit');
        if choice==0 || choice==4 || choice==5
            exportAnnotated();
            if isgraphics(fig), close(fig); end
            if choice==5, fprintf('All measurements saved.\n'); return; end
            break;
        elseif choice==1
            answer=inputdlg('Target point ID:','Add image point',[1 50]);
            if isempty(answer) || isempty(strtrim(answer{1})),continue;end
            id=string(strtrim(answer{1}));
            if contains(id,whitespacePattern)
                warndlg('IDs cannot contain whitespace.');continue;
            end
            if any(T.image_id==idImage & T.point_id==id)
                warndlg('That point is already measured in this image. Use Edit.');continue;
            end
            pos=selectPoint(id,[]);
            if isempty(pos),continue;end
            T(end+1,:)={id,idImage,imagePath,pos(1),pos(2)};
            saveData();refreshLabels();
        elseif choice==2 || choice==3
            idx=find(T.image_id==idImage);
            if isempty(idx)
                warndlg('No points measured in this image.');continue;
            end
            labels=cellstr(T.point_id(idx));
            [selected,ok]=listdlg('PromptString','Select point ID:', ...
                'SelectionMode','single','ListString',labels);
            if ~ok,continue;end
            j=idx(selected);
            if choice==2
                pos=selectPoint(T.point_id(j),[T.x_pix(j),T.y_pix(j)]);
                if isempty(pos),continue;end
                T.x_pix(j)=pos(1);T.y_pix(j)=pos(2);
            else
                answer=questdlg(sprintf('Delete point %s from image %d?', ...
                    T.point_id(j),idImage),'Confirm delete','Delete','Cancel','Cancel');
                if ~strcmp(answer,'Delete'),continue;end
                T(j,:)=[];
            end
            saveData();refreshLabels();
        end
    end
end
fprintf('Finished. All measurements saved to %s\n',projectFolder);

    function pos=selectPoint(id,oldPos)
        pos=[];
        if ~isgraphics(fig),return;end
        title(ax,sprintf('Point %s: zoom/pan, then press PLACE TARGET',id), ...
            'Interpreter','none');
        zoom(fig,'on');
        btn=uicontrol(fig,'Style','pushbutton','String','PLACE TARGET', ...
            'Units','normalized','Position',[0.73 0.015 0.25 0.065], ...
            'FontSize',11,'Callback',@(~,~) uiresume(fig));
        uiwait(fig);
        if isgraphics(btn),delete(btn);end
        if ~isgraphics(fig),return;end
        zoom(fig,'off');pan(fig,'off');
        title(ax,sprintf('Point %s: drag marker, double-click to accept',id), ...
            'Interpreter','none');
        try
            if isempty(oldPos)
                p=drawpoint(ax,'Color','y','Label',char(id));
            else
                p=drawpoint(ax,'Position',oldPos,'Color','y','Label',char(id));
            end
            wait(p);
            if isvalid(p)
                testPos=p.Position;
                delete(p);
            else
                return;
            end
        catch ME
            warning('Selection cancelled/failed: %s',ME.message);
            return;
        end
        if numel(testPos)~=2 || any(~isfinite(testPos)) || ...
                testPos(1)<1 || testPos(2)<1 || ...
                testPos(1)>size(im,2) || testPos(2)>size(im,1)
            warndlg('Selected point lies outside the image.');return;
        end
        pos=testPos;
    end

    function refreshLabels()
        if ~isgraphics(ax),return;end
        delete(findobj(ax,'Tag','measuredPoint'));
        delete(findobj(ax,'Tag','measuredLabel'));
        idx=find(T.image_id==idImage);
        for n=reshape(idx,1,[])
            plot(ax,T.x_pix(n),T.y_pix(n),'r+', ...
                'MarkerSize',9,'LineWidth',1.2,'Tag','measuredPoint');
            text(ax,T.x_pix(n)+9,T.y_pix(n),char(T.point_id(n)), ...
                'Color','yellow','FontSize',9,'FontWeight','bold', ...
                'Interpreter','none','Tag','measuredLabel');
        end
        title(ax,sprintf('Image %d — %d recorded targets',idImage,numel(idx)));
    end

    function exportAnnotated()
        if ~isgraphics(ax),return;end
        % Temporarily show the full photograph for annotation exports.
        oldXLim = xlim(ax); oldYLim = ylim(ax);
        oldXMode = ax.XLimMode; oldYMode = ax.YLimMode;
        restoreView = onCleanup(@() restoreLimits(oldXLim,oldYLim,oldXMode,oldYMode));
        xlim(ax,[0.5 size(im,2)+0.5]);
        ylim(ax,[0.5 size(im,1)+0.5]);
        try
            exportgraphics(ax,fullfile(projectFolder, ...
                sprintf('image_%02d_annotated.png',idImage)), ...
                'Resolution',150);
        catch ME
            warning('Annotated image export failed: %s',ME.message);
        end
        clear restoreView;
    end

    function restoreLimits(oldXLim,oldYLim,oldXMode,oldYMode)
        if isgraphics(ax)
            xlim(ax,oldXLim); ylim(ax,oldYLim);
            ax.XLimMode = oldXMode; ax.YLimMode = oldYMode;
        end
    end

    function saveData()
        % Export CSV first, then bundle-ready file, and image mapping.
        writetable(T,pointsPath);
        writetable(L,lookupPath);
        pathPHO=fullfile(projectFolder,'calibration.PHO');
        [fid,msg]=fopen(pathPHO,'w');
        if fid<0,error('Unable to write PHO file: %s',msg);end
        c=onCleanup(@() fclose(fid));
        [~,order]=sortrows([T.image_id,(1:height(T))']);
        for ii=reshape(order,1,[])
            fprintf(fid,'%s %d %.2f %.2f\n',char(T.point_id(ii)), ...
                T.image_id(ii),T.x_pix(ii),T.y_pix(ii));
        end
        clear c;
        fprintf('Saved %d observations across %d registered images.\n', ...
            height(T),height(L));
    end
end
