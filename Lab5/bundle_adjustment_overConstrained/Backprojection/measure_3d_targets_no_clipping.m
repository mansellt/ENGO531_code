  function measure_3d_targets_no_clipping()
% MEASURE_3D_TARGETS  Pick labelled 3-D points from an XYZ point cloud.
%
% Input: TXT/CSV containing XYZ or XYZ RGB Intensity Nx Ny Nz.
%        RGB values may be 0..255 or 0..1. Original RGB is displayed.
% Output: target_coordinates.csv, columns point_id,X,Y,Z,cloud_index.
%         The CSV is updated immediately after every add/edit/delete.
%
% Requires a desktop MATLAB session with interactive figures.
% No Computer Vision or Lidar Toolbox required.
%
% Workflow: rotate, zoom with wheel/ZOOM +/-, and pan via arrow buttons.
% Do not use MATLAB toolbar zoom or pan; those change camera state.
% Do not enable MATLAB's built-in zoom tool; it may reintroduce clipping.
% When picking, click on a cloud point to create a MATLAB data tip,
% then click ACCEPT SELECTED POINT in the figure. The coordinates are
% retrieved from the ORIGINAL cloud point, not estimated from the screen.

    [name, folder] = uigetfile({'*.txt;*.csv','XYZ point clouds (*.txt, *.csv)'; ...
                              '*.*','All files'}, ...
                              'Select XYZ or XYZ+RGB point cloud');
    if isequal(name,0), return; end
    cloudPath = fullfile(folder,name);
    % Determine the separator from the first data line. The supplied
    % ENGO531 file uses spaces and a comment header; the previous cloud
    % uses commas. Both are supported without changing files.
    fid = fopen(cloudPath,'r');
    if fid<0, error('Cannot open point cloud: %s',cloudPath); end
    firstLine = fgetl(fid);
    secondLine = fgetl(fid);
    fclose(fid);
    if ~ischar(firstLine) || ~ischar(secondLine)
        error('The point-cloud file has no data rows.');
    end
    if contains(secondLine,',')
        delimiter = ',';
    else
        delimiter = ' ';
    end
    data = readmatrix(cloudPath,'FileType','text', ...
                      'Delimiter',delimiter,'NumHeaderLines',1);
    if size(data,2)<3
        error('Expected at least X, Y, Z in columns 1-3.');
    end
    % Preserve original data row numbers, including any removed invalid rows.
    originalIndices = (1:size(data,1))';
    valid = all(isfinite(data(:,1:3)),2);
    if ~all(valid)
        warning('Ignoring %d rows with invalid XYZ.',sum(~valid));
    end
    data = data(valid,:);
    originalIndices = originalIndices(valid);
    XYZ = data(:,1:3);
    if isempty(XYZ), error('No valid XYZ points were found.'); end

    hasRGB = size(data,2)>=6 && all(all(isfinite(data(:,4:6))));
    if hasRGB
        RGB = data(:,4:6);
        if any(RGB(:)>1)
            RGB = RGB/255;
        end
        RGB = min(1,max(0,RGB));
        fprintf('Detected %d columns; displaying original RGB from columns 4-6.\n',size(data,2));
    else
        RGB = [];
        fprintf('No valid RGB columns found; using a single cloud colour.\n');
    end
    fprintf('Loaded %d XYZ points from %s\n',size(XYZ,1),cloudPath);

    outputFolder = uigetdir(folder,'Choose folder for target coordinate CSV');
    if isequal(outputFolder,0), return; end
    csvPath = fullfile(outputFolder,'target_coordinates.csv');
    vars = {'point_id','X','Y','Z','cloud_index'};
    if isfile(csvPath)
        T = readtable(csvPath,'TextType','string');
        if ~all(ismember(vars,T.Properties.VariableNames))
            error('Existing target_coordinates.csv has unexpected columns. Back it up or choose another folder.');
        end
        T = T(:,vars);
        T.point_id = string(T.point_id);
        if any(~isfinite(T.cloud_index)) || any(T.cloud_index<1) || ...
                any(~ismember(T.cloud_index,originalIndices)) || any(T.cloud_index~=fix(T.cloud_index))
            error('Existing CSV has invalid cloud_index values for this cloud.');
        end
        if numel(unique(T.point_id)) ~= height(T)
            error('The CSV has duplicate point IDs. Resolve them before proceeding.');
        end
        [~,previousRows] = ismember(T.cloud_index,originalIndices);
        if any(any(abs([T.X,T.Y,T.Z]-XYZ(previousRows,:))>1e-7))
            error(['Existing CSV coordinates do not match this input cloud. ' ...
                   'Choose the original cloud or a different output folder.']);
        end
    else
        T = table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
                  zeros(0,1),'VariableNames',vars);
    end
    fprintf('Loaded %d previously recorded targets.\n',height(T));

    fig = figure('Name','3D Target Coordinate Measurement', ...
                 'NumberTitle','off','Color','w','Position',[80 80 1150 750]);
    ax = axes('Parent',fig,'Position',[0.07 0.12 0.86 0.80]);
    if hasRGB
        cloudColors = RGB;   % One true-colour [R G B] triplet per point
    else
        cloudColors = [0.32 0.48 0.58];
    end
    hCloud = scatter3(ax,XYZ(:,1),XYZ(:,2),XYZ(:,3), ...
                      3,cloudColors,'.');
    hCloud.PickableParts = 'visible';
    % Disable axes-boundary clipping, including on the point cloud itself.
    % Otherwise a close-up or translated view may show a sharp cut-off.
    ax.Clipping = 'off';
    hCloud.Clipping = 'off';
    hold(ax,'on');
    axis(ax,'equal'); axis(ax,'vis3d');
    grid(ax,'on'); xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
    view(ax,3);
    set(ax,'Projection','orthographic');

    % CUSTOM CLOSE-UP ZOOM: Keep the camera safely OUTSIDE the cloud.
    % MATLAB's normal scroll-wheel zoom can move its camera into the cloud,
    % causing near-plane clipping even with orthographic projection.
    % Instead, change only CameraViewAngle; CameraPosition stays fixed.
    xyzMin = min(XYZ,[],1);
    xyzMax = max(XYZ,[],1);
    cloudCentre = (xyzMin + xyzMax)/2;
    cloudRadius = max(norm(xyzMax-xyzMin)/2,1e-6);
    % Maintain fixed, generously padded data limits while navigating.
    % Turning clipping OFF is the critical part; the limits are for scale.
    margin = max(xyzMax-xyzMin,cloudRadius*0.05)*0.5;
    xlim(ax,[xyzMin(1)-margin(1),xyzMax(1)+margin(1)]);
    ylim(ax,[xyzMin(2)-margin(2),xyzMax(2)+margin(2)]);
    zlim(ax,[xyzMin(3)-margin(3),xyzMax(3)+margin(3)]);
    set(ax,'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
    viewingDirection = ax.CameraPosition - ax.CameraTarget;
    viewingDirection = viewingDirection/norm(viewingDirection);
    ax.CameraTarget = cloudCentre;
    ax.CameraPosition = cloudCentre + viewingDirection*(25*cloudRadius);
    initialViewAngle = 2*atan2d(1.45*cloudRadius,25*cloudRadius);
    ax.CameraViewAngle = initialViewAngle;
    ax.CameraViewAngleMode = 'manual';
    ax.CameraPositionMode = 'manual';
    ax.CameraTargetMode = 'manual';
    rotate3d(fig,'on');

    % Override wheel zoom; do NOT use MATLAB's standard zoom tool.
    fig.WindowScrollWheelFcn = @customWheelZoom;
    uicontrol(fig,'Style','pushbutton','String','ZOOM +', ...
        'Units','normalized','Position',[0.08 0.018 0.12 0.05], ...
        'Callback',@(~,~) changeZoom(-1));
    uicontrol(fig,'Style','pushbutton','String','ZOOM -', ...
        'Units','normalized','Position',[0.21 0.018 0.12 0.05], ...
        'Callback',@(~,~) changeZoom(1));
    % Pan/translate by moving the camera and target together; this keeps
    % the camera far outside the cloud (unlike built-in camera translation).
    uicontrol(fig,'Style','pushbutton','String','LEFT', ...
        'Units','normalized','Position',[0.49 0.018 0.08 0.05], ...
        'Callback',@(~,~) panView(-1,0));
    uicontrol(fig,'Style','pushbutton','String','RIGHT', ...
        'Units','normalized','Position',[0.58 0.018 0.08 0.05], ...
        'Callback',@(~,~) panView(1,0));
    uicontrol(fig,'Style','pushbutton','String','UP', ...
        'Units','normalized','Position',[0.67 0.018 0.07 0.05], ...
        'Callback',@(~,~) panView(0,1));
    uicontrol(fig,'Style','pushbutton','String','DOWN', ...
        'Units','normalized','Position',[0.75 0.018 0.08 0.05], ...
        'Callback',@(~,~) panView(0,-1));
    uicontrol(fig,'Style','pushbutton','String','RESET VIEW', ...
        'Units','normalized','Position',[0.84 0.018 0.13 0.05], ...
        'Callback',@(~,~) resetZoom());
    refreshLabels();
    % MATLAB may reenable clipping after certain graphics interactions.
    ax.Clipping = 'off';
    hCloud.Clipping = 'off';

    while isgraphics(fig)
        choice = menu(sprintf('%d targets recorded | Use Rotate/zoom view before selecting',height(T)), ...
            'Rotate / zoom view','Add target','Edit existing target', ...
            'Delete target','Save and close');
        if choice==0 || choice==5
            if isgraphics(fig), close(fig); end
            break;
        end
        if choice==1
            title(ax,'Rotate/zoom the point cloud, then press DONE VIEWING');
            doneButton = uicontrol(fig,'Style','pushbutton', ...
                'String','DONE VIEWING','Units','normalized', ...
                'Position',[0.73 0.018 0.24 0.055], ...
                'Callback',@(~,~) uiresume(fig));
            uiwait(fig);
            if isgraphics(doneButton), delete(doneButton); end
            if isgraphics(fig), refreshLabels(); end
            continue;
        end
        if choice==2
            answer = inputdlg('Enter unique target ID (e.g., P1):', ...
                              'New 3D target',[1 50]);
            if isempty(answer), continue; end
            id = strtrim(string(answer{1}));
            if strlength(id)==0, continue; end
            if any(T.point_id==id)
                warndlg('That target already exists. Choose Edit instead.');
                continue;
            end
            pointIndex = pickCloudPoint(id);
            if isempty(pointIndex), continue; end
            p = XYZ(pointIndex,:);
            T(end+1,:) = {id,p(1),p(2),p(3),originalIndices(pointIndex)};
            saveCSV(); refreshLabels();
        elseif choice==3 || choice==4
            if isempty(T)
                warndlg('No points have been recorded yet.'); continue;
            end
            [idx,ok] = listdlg('ListString',cellstr(T.point_id), ...
                'SelectionMode','single','PromptString','Select target ID:', ...
                'ListSize',[300 320]);
            if ~ok, continue; end
            if choice==3
                pointIndex = pickCloudPoint(T.point_id(idx));
                if isempty(pointIndex), continue; end
                p = XYZ(pointIndex,:);
                T.X(idx)=p(1); T.Y(idx)=p(2); T.Z(idx)=p(3);
                T.cloud_index(idx)=originalIndices(pointIndex);
            else
                answer = questdlg(sprintf('Delete target %s?',T.point_id(idx)), ...
                    'Confirm deletion','Delete','Cancel','Cancel');
                if ~strcmp(answer,'Delete'), continue; end
                T(idx,:)=[];
            end
            saveCSV(); refreshLabels();
        end
    end
    fprintf('Finished. CSV location: %s\n',csvPath);

    function customWheelZoom(~,evt)
        % VerticalScrollCount < 0 means zoom in on most mouse wheels.
        if ~isgraphics(ax), return; end
        changeZoom(evt.VerticalScrollCount);
    end

    function changeZoom(steps)
        if ~isgraphics(ax), return; end
        % Only the viewing angle changes, never the camera distance.
        ax.CameraViewAngle = min(90,max(0.005, ...
            ax.CameraViewAngle*(1.25^steps)));
        drawnow limitrate;
    end

    function panView(horizontal,vertical)
        if ~isgraphics(ax), return; end
        forward = ax.CameraTarget - ax.CameraPosition;
        forward = forward/norm(forward);
        up = ax.CameraUpVector;
        up = up/norm(up);
        right = cross(forward,up);
        right = right/norm(right);
        up = cross(right,forward);
        % Shift 15 percent of the currently visible half-height.
        visibleHalfHeight = norm(ax.CameraPosition-ax.CameraTarget) * ...
                            tand(ax.CameraViewAngle/2);
        delta = 0.3*visibleHalfHeight * ...
                (-horizontal*right - vertical*up);
        ax.CameraTarget = ax.CameraTarget + delta;
        ax.CameraPosition = ax.CameraPosition + delta;
        drawnow limitrate;
    end

    function resetZoom()
        if ~isgraphics(ax), return; end
        ax.CameraViewAngle = initialViewAngle;
        ax.CameraTarget = cloudCentre;
        ax.CameraPosition = cloudCentre + viewingDirection*(25*cloudRadius);
        drawnow;
    end

    function pointIndex = pickCloudPoint(id)
        % A data-tip callback retains the exact selected source point.
        % This avoids relying on getCursorInfo when ACCEPT is pressed.
        pointIndex = [];
        if ~isgraphics(fig), return; end
        rotate3d(fig,'off');
        dcm = datacursormode(fig);
        set(dcm,'Enable','off');
        removeAllDataTips(dcm);
        selectedIndex = [];
        set(dcm,'UpdateFcn',@onDataTip, ...
                'DisplayStyle','datatip','SnapToDataVertex','on', ...
                'Enable','on');
        title(ax,sprintf('CLICK one cloud point for %s, then ACCEPT SELECTED POINT',id), ...
              'Interpreter','none');
        accept = uicontrol(fig,'Style','pushbutton', ...
            'String','ACCEPT SELECTED POINT', ...
            'Units','normalized','Position',[0.70 0.018 0.27 0.055], ...
            'FontSize',11,'BackgroundColor',[0.86 0.95 0.85], ...
            'Callback',@(~,~) acceptPick());
        cancel = uicontrol(fig,'Style','pushbutton','String','CANCEL', ...
            'Units','normalized','Position',[0.51 0.018 0.17 0.055], ...
            'Callback',@(~,~) cancelPick());
        uiwait(fig);
        if isgraphics(accept), delete(accept); end
        if isgraphics(cancel), delete(cancel); end
        if isgraphics(fig)
            set(dcm,'Enable','off');
            removeAllDataTips(dcm);
            rotate3d(fig,'on');
            refreshLabels();
        end

        function output = onDataTip(~,eventObj)
            pos = get(eventObj,'Position');
            try
                dataIdx = get(eventObj,'DataIndex');
            catch
                dataIdx = [];
            end
            if isempty(dataIdx) || ~isscalar(dataIdx) || ...
                    dataIdx<1 || dataIdx>size(XYZ,1)
                % Fall back to matching the snapped XYZ coordinate.
                [~,dataIdx] = min(sum((XYZ-pos).^2,2));
            end
            selectedIndex = dataIdx;
            q = XYZ(dataIdx,:);
            output = {sprintf('Target: %s',char(id)), ...
                      sprintf('Cloud index: %d',originalIndices(dataIdx)), ...
                      sprintf('X: %.9f',q(1)), ...
                      sprintf('Y: %.9f',q(2)), ...
                      sprintf('Z: %.9f',q(3))};
            fprintf('Selected cloud point %d for target %s\n', ...
                    originalIndices(dataIdx),char(id));
        end

        function acceptPick()
            if isempty(selectedIndex)
                % Some MATLAB releases update data tips differently.
                % Try the active data tip before asking for another click.
                try
                    tips = getCursorInfo(dcm);
                    if ~isempty(tips)
                        matching = find(arrayfun(@(x) isequal(x.Target,hCloud),tips),1,'end');
                        if ~isempty(matching)
                            selectedInfo = tips(matching);
                            if isfield(selectedInfo,'DataIndex') && ...
                                    ~isempty(selectedInfo.DataIndex)
                                selectedIndex = selectedInfo.DataIndex;
                            else
                                pos = selectedInfo.Position;
                                [~,selectedIndex] = min(sum((XYZ-pos).^2,2));
                            end
                        end
                    end
                catch
                end
            end
            if isempty(selectedIndex)
                warndlg(['No point selected yet. Click directly on a dot ' ...
                         'in the cloud until an XYZ data tip appears, then press ACCEPT.'], ...
                         'Select a cloud point');
                return;
            end
            q = XYZ(selectedIndex,:);
            answer = questdlg(sprintf(['Target: %s\n' ...
                  'X = %.9f\nY = %.9f\nZ = %.9f\nCloud point = %d\n\nRecord these coordinates?'], ...
                   char(id),q(1),q(2),q(3),originalIndices(selectedIndex)), ...
                   'Confirm 3D coordinates','Record','Pick again','Cancel','Record');
            if strcmp(answer,'Record')
                pointIndex = selectedIndex;
                uiresume(fig);
            elseif strcmp(answer,'Cancel')
                pointIndex = [];
                uiresume(fig);
            else
                selectedIndex = [];
                removeAllDataTips(dcm);
            end
        end

        function cancelPick()
            pointIndex = [];
            if isgraphics(fig), uiresume(fig); end
        end
    end

    function refreshLabels()
        if ~isgraphics(ax), return; end
        delete(findobj(ax,'Tag','recordedTarget'));
        delete(findobj(ax,'Tag','recordedLabel'));
        for j=1:height(T)
            plot3(ax,T.X(j),T.Y(j),T.Z(j),'ro','MarkerSize',9, ...
                  'LineWidth',1.8,'MarkerFaceColor','y','Tag','recordedTarget');
            text(ax,T.X(j),T.Y(j),T.Z(j),['  ' char(T.point_id(j))], ...
                 'Color','r','FontSize',11,'FontWeight','bold', ...
                 'Interpreter','none','Tag','recordedLabel');
        end
        title(ax,sprintf('%d cloud points | %d recorded targets',size(XYZ,1),height(T)));
        ax.Clipping = 'off';
        hCloud.Clipping = 'off';
    end

    function saveCSV()
        % Temp file then replacement: helps avoid a partly written CSV.
        tmpPath = [csvPath '.tmp'];
        writetable(T,tmpPath,'FileType','text');
        [ok,msg] = movefile(tmpPath,csvPath,'f');
        if ~ok, error('Could not save CSV: %s',msg); end
        fprintf('Saved %d targets to %s\n',height(T),csvPath);
    end
end

function removeAllDataTips(dcm)
% Delete only data cursor objects, not the original scatter plot.
try
    tips = getCursorInfo(dcm);
    for j=1:numel(tips)
        if isfield(tips(j),'DataCursor') && isgraphics(tips(j).DataCursor)
            delete(tips(j).DataCursor);
        end
    end
catch
    % Some MATLAB releases manage data-tip deletion automatically.
end
end
