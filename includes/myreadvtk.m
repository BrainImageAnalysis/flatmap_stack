function datasets = myreadvtk(fname)
%

fid = fopen(fname);
tline = fgetl(fid);
state = "DATASET POLYDATA";
state = "";
substate = 0;
datasets = [];
while ischar(tline)
  %  disp(tline)
   switch (state)
       case "DATASET_POLYDATA"
            tline_data = strsplit(tline,' ');
            switch (substate)
                case 0
                    fprintf('substate DATASET POLYDATA\n');
                    num_pts = str2num(tline_data{2});
                    points = zeros(num_pts*3,1);
                    V = [];                
                    datasets = setfield(datasets,state,V); 

                    fprintf('%d %d\n',size(points));
                    substate = 1;
                    datasets.DATASET_POLYDATA.points = points;
                    points_offset = 1;
                case 1
                    if ~strcmp(tline_data{1},'POLYGONS')
                        points = str2num(tline);
                        npoints = numel(points);
                        datasets.DATASET_POLYDATA.points(points_offset:points_offset+npoints-1) = points;
                        points_offset = points_offset + npoints;
                    else
                         datasets.DATASET_POLYDATA.points = reshape(datasets.DATASET_POLYDATA.points,[3,numel(datasets.DATASET_POLYDATA.points)/3]);
                         fprintf('%s\n',tline);
                         num_pts = str2num(tline_data{2});
                         points = zeros(num_pts*3,1);
                         fprintf('%d %d\n',size(points));
                         V = [];   
                         substate = 2;
                         datasets.DATASET_POLYDATA.polygons = points;
                         points_offset = 1;
                    end
                case 2
                    if ~strcmp(tline_data{1},'POINT_DATA')
                        points = str2num(tline);
                        points = points(2:end);
                        npoints = numel(points);
                        datasets.DATASET_POLYDATA.polygons(points_offset:points_offset+npoints-1) = points;
                        points_offset = points_offset + npoints;
                    else
                         datasets.DATASET_POLYDATA.polygons = reshape(datasets.DATASET_POLYDATA.polygons,[3,numel(datasets.DATASET_POLYDATA.polygons)/3]);
                         V = [];   
                         substate = 3;
                    end
                case 3
                    if strcmp(tline_data{1},'NORMALS')
                        substate = 4;             
                      %  fprintf('%s\n',tline);
                         num_pts = numel(datasets.DATASET_POLYDATA.points);
                         points = zeros(num_pts,1);
                         fprintf('%d %d\n',size(points));
                         V = [];   
                         %substate = 2;
                         datasets.DATASET_POLYDATA.normals = points;
                         points_offset = 1;
                    end
                case 4
                    if ~strcmp(tline_data{1},'FIELD')
                        %if ~strcmp(tline_data{1},'Label_ID')
                        points = str2num(tline);
                        npoints = numel(points);
                        datasets.DATASET_POLYDATA.normals(points_offset:points_offset+npoints-1) = points;
                        points_offset = points_offset + npoints;                     
                        %end
                    else
                         datasets.DATASET_POLYDATA.normals = reshape(datasets.DATASET_POLYDATA.normals,[3,numel(datasets.DATASET_POLYDATA.normals)/3]);
                         fprintf('%s\n',tline);
                         substate = 6;
                         
                         %substate = 6;             
                         num_pts = numel(datasets.DATASET_POLYDATA.points);
                         points = zeros(num_pts/3,1);
                         fprintf('%d %d\n',size(points));
                         V = [];   
                         datasets.DATASET_POLYDATA.atlas = points;
                         points_offset = 1;
                         
                    end
                    
                
                case 6
                    if ~strcmp(tline_data{1},'unknown')
                        if ~strcmp(tline_data{1},'Label_ID')
                        points = str2num(tline);
                        npoints = numel(points);
                        datasets.DATASET_POLYDATA.atlas(points_offset:points_offset+npoints-1) = points;
                        points_offset = points_offset + npoints;
                        end
                    else
                        %datasets.DATASET_POLYDATA.atlas = reshape(datasets.DATASET_POLYDATA.atlas,[3,numel(datasets.DATASET_POLYDATA.atlas)/3]);
                        substate = 7; 
                    end
               
                   
            end
            
            
   end
    switch (tline)
        case "DATASET POLYDATA"
            state = "DATASET_POLYDATA";
            fprintf('DATASET POLYDATA\n');
            substate = 0;
    
            
    end
    tline = fgetl(fid);
end
fclose(fid);

%%
if false
    %%
    datasets = myreadvtk('brain_flatmap.vtk');
    datasets_3D = myreadvtk('brain_midthickness.vtk');

    %%
    gs = [];
    clf
    gs.faces =  datasets.DATASET_POLYDATA.polygons'+1;
    gs.vertices = datasets.DATASET_POLYDATA.points';

    colors = squeeze(cat(2,datasets.DATASET_POLYDATA.atlas,datasets.DATASET_POLYDATA.atlas,datasets.DATASET_POLYDATA.atlas));

    %patch(gs,'FaceColor','blue','EdgeColor','none');FaceVertexCData
    %colors = rand(size(colors));
    patch(gs,'FaceVertexCData',colors/max(colors(:)),'FaceColor','interp','EdgeColor','none');
    daspect([1,1,.4]); view(45,30); axis tight

    lightangle(45,30);

end


