
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   global settings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
highres = false;
%FM = 'BM';

FM = 'MBM';
%highres = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label2 = load_untouch_nii('data/cortex_boundary_TC_std.nii.gz');
label2 = label2.img;
load('data/fmstack_data_n.mat','tracer_mask','tracer_gradf','cortex_dist');
%load('~/tmp/tracer_data_n.mat','tracer_mask','tracer_gradf','cortex_dist');
%%


out_postfix = '';


%%

TC_img = load_untouch_nii(TC_avg_nii_fn);    
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   flatstack tracing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




min_sampling_pts = 3;

switch FM
    
    case 'BM'
        ofolder = [debug_out_fd,'/BMA/'];
        %gi_surf = gifti('data/surf/BMA/midsurface_2018_left_org.surf.gii');
        gi_surf = gifti('data/surf/BMA/bma_sp2.lh.graymid.surf.gii');
        gi_flat = gifti('data/surf/BMA/bma_sp2.lh.flatmap.surf.gii');
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   reference img (Tc_ standard) for mapping the coordinates
        %   ti TC_std brain voxel space
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        addpath data/surf/BMA/
   
        if false
           verts_tc_std = map_pts_BMA(gi_surf.vertices')';
           save('data/surf/BMA/pts_TC.mat','verts_tc_std');
        else
           %use precomputed points 
           load('data/surf/BMA/pts_TC.mat','verts_tc_std');
        end
        gi_surf2 = gi_surf;
        gi_surf2.vertices(:) = verts_tc_std(:);
    
        %%
        hdr = TC_img.hdr;
        T = [hdr.hist.srow_x;hdr.hist.srow_y;hdr.hist.srow_z;0,0,0,1];
        T = pinv(T);
        pts = [gi_surf2.vertices,ones(size(gi_surf2.vertices,1),1)];
        pts = (T*pts')';
        gi_surf_tc_vox = gi_surf2;
        gi_surf_tc_vox.vertices(:) = pts(:,1:3);
        shape = size(TC_img.img);
        TC_img_midsurf = zeros(shape,'single');
    case 'MBM'
        ofolder = [debug_out_fd,'/MBM/'];
        mbm = load('data/surf/MBM/MBM_fmap.mat');
        gi_surf = mbm.surf{2};
        gi_flat = mbm.fm;    
        pts = gi_surf.vertices;
        
        gi_surf_tc_vox = gi_surf;
        shape = size(TC_img.img);
        TC_img_midsurf = zeros(shape,'single');

        if highres 
            ofolder = [debug_out_fd,'/MBM_highres/'];
            min_sampling_pts = 2;
            dist_check = 1;
            safe_margin = 0.1;
            
            dist_check = 5;
            safe_margin = 0.05;
            safe_margin = 0.1;
            
            old_verts1 = gi_surf.vertices(gi_surf.faces(:,1),:);
            old_verts2 = gi_surf.vertices(gi_surf.faces(:,2),:);
            old_verts3 = gi_surf.vertices(gi_surf.faces(:,3),:);

            dist = @(a,b)sqrt(sum((a-b).^2,2));
            new_verts = (old_verts1+old_verts2+old_verts3)/3.0;
            vert_dist = max(max(dist(old_verts2,old_verts1),dist(old_verts1,old_verts3)),dist(old_verts3,old_verts2));


            old_verts1_r = gi_flat.vertices(gi_surf.faces(:,1),:);
            old_verts2_r = gi_flat.vertices(gi_surf.faces(:,2),:);
            old_verts3_r = gi_flat.vertices(gi_surf.faces(:,3),:);
            new_verts_r = (old_verts1+old_verts2+old_verts3)/3.0;
            vert_dist2 = max(max(dist(old_verts2_r,old_verts1_r),dist(old_verts1_r,old_verts3_r)),dist(old_verts3_r,old_verts2_r));
            flatmap_safe = vert_dist2<dist_check;


            
            safezone = ((0.5-abs(cortex_dist-0.5))>safe_margin*0.5).*(label2==1);
            safezone_dist = (vert_dist<5) & flatmap_safe;

            new_verts_2 = (new_verts+old_verts2+old_verts3)/3.0;
            new_verts_3 = (old_verts1+new_verts+old_verts3)/3.0;
            new_verts_4 = (old_verts1+old_verts2+new_verts)/3.0;

            new_verts = cat(1,new_verts(safezone_dist,:),new_verts_2(safezone_dist,:),new_verts_3(safezone_dist,:),new_verts_4(safezone_dist,:));

            indx = sub2ind(shape,round(new_verts(:,1)),round(new_verts(:,2)),round(new_verts(:,3)));
            new_valid_verts = (safezone(indx)>0);

            indx2 = sub2ind(shape,round(gi_surf_tc_vox.vertices(:,1)),round(gi_surf_tc_vox.vertices(:,2)),round(gi_surf_tc_vox.vertices(:,3)));
            old_valid_verts = (safezone(indx2)>0);

            gi_surf_tc_vox.vertices = cat(1,new_verts(new_valid_verts,:),gi_surf_tc_vox.vertices(old_valid_verts,:));

            gi_surf.vertices = gi_surf_tc_vox.vertices;
            pts = gi_surf.vertices;
            if true
                old_verts1 = gi_flat.vertices(gi_surf.faces(:,1),:);
                old_verts2 = gi_flat.vertices(gi_surf.faces(:,2),:);
                old_verts3 = gi_flat.vertices(gi_surf.faces(:,3),:);

                new_verts = (old_verts1+old_verts2+old_verts3)/3.0;
                new_verts_2 = (new_verts+old_verts2+old_verts3)/3.0;
                new_verts_3 = (old_verts1+new_verts+old_verts3)/3.0;
                new_verts_4 = (old_verts1+old_verts2+new_verts)/3.0;            

                new_verts = cat(1,new_verts(safezone_dist,:),new_verts_2(safezone_dist,:),new_verts_3(safezone_dist,:),new_verts_4(safezone_dist,:));
                gi_flat.vertices = cat(1,new_verts(new_valid_verts,:),gi_flat.vertices(old_valid_verts,:));
            end

        end
end

mkdir(ofolder);
indx = sub2ind(shape,round(pts(:,1)),round(pts(:,2)),round(pts(:,3)));

TC_img_midsurf(:) = 0;
TC_img_midsurf(indx)=1;

trace_nii = TC_img;
trace_nii.hdr.dime.datatype = 128;
trace_nii.hdr.dime.bitpix = 24;
trace_nii.hdr.dime.glmax = 255;
trace_nii.hdr.dime.glmin = 0;
trace_nii.img = uint8(cat(4,label2,TC_img_midsurf,TC_img_midsurf)*85);

save_untouch_nii(trace_nii,[ofolder,'/debug01_center_pts.nii.gz']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   trace the cortex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[result_tracer,result_tracks] = trace_flatmap_stack({'ofield',tracer_gradf,'mask',tracer_mask,'seeds',pts(:,[1:3]),'dist',cortex_dist,'density',10});


%%
 
trace_nii = TC_img;
trace_nii.hdr.dime.datatype = 128;
trace_nii.hdr.dime.bitpix = 24;
trace_nii.hdr.dime.glmax = 255;
trace_nii.hdr.dime.glmin = 0;
trace_nii.img = uint8(permute(result_tracer*255,[2,3,4,1]));

save_untouch_nii(trace_nii,[ofolder,'/debug02_traces.nii.gz']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   resample traces equidistantly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mode = 'equidist';
%mode = 'heat';


smpl_pts = 50;
intp = [0:smpl_pts-1]/smpl_pts;
result_tracks_array = zeros([numel(result_tracks),smpl_pts,3]);
for a = 1:numel(result_tracks)
    switch mode 
        case 'heat'
            data = result_tracks{a};
            if (size(data,1)>3)
                [uval,uind] = unique(squeeze(data(:,4)));
                
                if numel(uind)>3
                uval(1) = 0;
                uval(end) = 1;
                vq = interp1(uval,data(uind,1:3),intp);
                result_tracks_array(a,:,:) = vq;
                if any(isnan(vq))
                    disp(vq)
                    break
                end
                end
            end
        case 'equidist'

            data = result_tracks{a}(:,1:3);
            if (size(data,1)>min_sampling_pts)
            d = [0:size(data,1)-1]/(size(data,1)-1);

            vq = interp1(d,data,intp);
            result_tracks_array(a,:,:) = vq;
            end
    end
end


%%
if true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   map layers to TC_std img
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TC_img_midsurf = zeros(shape,'single');
    
    count = 1;
    valid = max(squeeze(result_tracks_array(:,:,1)),[],2) > 0;
    DI = (result_tracks_array(valid,2:end,:)-result_tracks_array(valid,1:end-1,:));
    DI = cat(2,zeros([size(DI,1),1]),sqrt(sum(DI.^2,3)));
    DIc = cumsum(DI,2);
    for di = [1,5,10,15,20,25,30,35,40,45,50]
        pts_l = squeeze(result_tracks_array(valid,di,:))+1;
        indx = sub2ind(shape,round(pts_l(:,1)),round(pts_l(:,2)),round(pts_l(:,3)));

        TC_img_midsurf(indx) = count;
       
        count = count + 1;
    end
    n = max(TC_img_midsurf(:));
    
    trace_nii = TC_img;
    trace_nii.hdr.dime.datatype = 128;
    trace_nii.hdr.dime.bitpix = 24;
    trace_nii.hdr.dime.glmax = 255;
    trace_nii.hdr.dime.glmin = 0;
    trace_nii.img =uint8(cat(4,label2*0.25*85,255*TC_img_midsurf/n,255*(TC_img_midsurf>0)));
    save_untouch_nii(trace_nii,[ofolder,'/debug03_equidistant_trace_samples.nii.gz']);
    
    save([ofolder,'dist.mat'],'DIc');
    %%
    
    iters = 5000;

    H = zeros([3,3,3],'double');
    H(2,2,1) = 1;
    H(:,:,2) = [0,1,0;1,-6,1;0,1,0];
    H(2,2,3) = 1;

    data = TC_img_midsurf;
    data(label2==3) = max(TC_img_midsurf(:));
    data(label2==2) = 1;

    mask = data>0;
    data2_org = gpuArray(data);


    data2 = gpuArray(data);
    mask2 = gpuArray(mask);
    H2 =  gpuArray(H);
    for a = 1:iters
        fprintf('%d \n',a);
        data2_lap = imfilter(data2,H2);
        data2 = data2 + 0.05 * data2_lap;
        data2(mask2) = data2_org(mask2);
        if mod(a-1,5) == 0
            sfigure(1);
            imagesc(squeeze(data2(:,:,ceil(end/2))));
            title(num2str(a));
            drawnow
        end
    end
    layer_img = gather(data2);
    %%
    trace_nii = TC_img;
    trace_nii.hdr.dime.datatype = 16;
    trace_nii.hdr.dime.bitpix = 32;
    trace_nii.img = single(layer_img);
    save_untouch_nii(trace_nii,[ofolder,'/debug04_equidistant_layers_continous.nii.gz']);

    
    
    %%

    layer_img_ = min(round((layer_img+0.5).*(label2>0)), max(TC_img_midsurf(:)));
   
    layer_img_ = medfilt3(layer_img_);

    uids = unique(layer_img_(:));
    uids(uids==0) = [];
    for a = 1:numel(uids)
    layer_img_(layer_img_==uids(a)) = a;
    end
 
    %%
    trace_nii = TC_img;
    trace_nii.hdr.dime.datatype = 2;
    trace_nii.hdr.dime.bitpix = 8;
    trace_nii.hdr.dime.glmax = 255;
    trace_nii.hdr.dime.cal_max = max(layer_img_(:));
    
    trace_nii.img = uint8(layer_img_);
    trace_nii.img(end:-1:end/2+1,:,:) = uint8(layer_img_(1:end/2,:,:));
    save_untouch_nii(trace_nii,[ofolder,'/debug05_layers_discrete.nii.gz']);
    %%
  
end
%%


if ~highres
    surfaces = {};
    meshes = {};

    for di = 1:size(result_tracks_array,2)
        fprintf('%d %d\n',di,size(result_tracks_array,2));
        valid_verts = max(squeeze(result_tracks_array(:,:,1)),[],2) > 0;

        pts_l = squeeze(result_tracks_array(valid_verts,di,:))+1;
        indx = sub2ind(shape,round(pts_l(:,1)),round(pts_l(:,2)),round(pts_l(:,3)));
        g_x = squeeze(tracer_gradf(1,:,:,:));
        g_y = squeeze(tracer_gradf(2,:,:,:));
        g_z = squeeze(tracer_gradf(3,:,:,:));
        surface_n = cat(2,g_x(indx),g_y(indx),g_z(indx));
        surface_n = surface_n./repmat(sqrt(sum(surface_n.^2,2)),1,3);

        valid_faces = valid_verts(gi_surf.faces(:,1)) & valid_verts(gi_surf.faces(:,2)) & valid_verts(gi_surf.faces(:,3));
        new_verts_indx =zeros(size(valid_verts));
        new_verts_indx(valid_verts) = 1:sum(valid_verts);
        gi_surf_new = gi_surf;
        gi_surf_new.vertices=pts_l;%gi_surf.vertices(valid_verts,:);
        for a = 1:3
            gi_surf_new.faces(valid_faces,a) = new_verts_indx(gi_surf.faces(valid_faces,a));
        end
        gi_surf_new.faces = gi_surf_new.faces(valid_faces,:);
        surfaces{di} = gi_surf_new;
        surfaces{di}.normals = surface_n;
        mesh = [];
        mesh.v = single(surfaces{di}.vertices)'-1;
        mesh.f = uint32(surfaces{di}.faces'-1);
        mesh.n = single(surfaces{di}.normals).';
        mesh.ni = mesh.f;
        meshes{di} = mesh;
    end

    switch FM
        case 'BM'
            save('./data/surfaces_BMA.mat','surfaces','meshes');
        case 'MBM'
            save('./data/surfaces_MBM.mat','surfaces','meshes');
    end
    
    %mesh_img = render_mesh2img({'faces',single(surfaces{mind}.faces-1),'vertices',single(surfaces{mind}.vertices),'shape',size(TC_img.img)});
    %start_VIannotator(mesh_img)
    %start_VIannotator({label2*0.25,mesh_img,TC_img_midsurf>0})
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  compute deformation (area/mid area ratio)
%  [computation continues some paragraphs below]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
valid = max(squeeze(result_tracks_array(:,:,1)),[],2) > 0;
faces = gi_surf_tc_vox.faces;
v1 = result_tracks_array(faces(:,1),:,:);
v2 = result_tracks_array(faces(:,2),:,:);
v3 = result_tracks_array(faces(:,3),:,:);
area = sqrt(sum(cross((v1-v2),(v3-v2)).^2,3))/2;
%%
vertex_feat = zeros(size(result_tracks_array,1),size(result_tracks_array,2));
for a = 1:size(faces,1)
    indx1 = faces(a,1);
    vertex_feat(indx1,:) = vertex_feat(indx1,:) + area(a,:);
    indx1 = faces(a,2);
    vertex_feat(indx1,:) = vertex_feat(indx1,:) + area(a,:);
    indx1 = faces(a,3);
    vertex_feat(indx1,:) = vertex_feat(indx1,:) + area(a,:);
end

vertex_feat = vertex_feat ./ repmat(vertex_feat(:,end/2),1,size(vertex_feat,2));
vertex_feat(isnan(vertex_feat)) = 0;
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  flatmap mapping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valid = max(squeeze(max(result_tracks_array(:,:,:),[],3)),[],2) > 0;

indx_3D = sub2ind(shape,round(pts(:,1)),round(pts(:,2)),round(pts(:,3)));
TC_img_midsurf(:) = 0;
TC_img_midsurf(indx_3D(valid))=0;
TC_img_midsurf(indx_3D(~valid))=1;

trace_nii = TC_img;
trace_nii.hdr.dime.datatype = 2;
trace_nii.hdr.dime.bitpix = 8;
trace_nii.hdr.dime.glmax = 255;
trace_nii.hdr.dime.glmin = 0;
trace_nii.img = uint8(TC_img_midsurf*255);
save_untouch_nii(trace_nii,[ofolder,'/debug06_missing.nii.gz']);
%%

f_shape = [1000,1000];


switch FM
    case 'Rosa'
        border = 0;
        offset = [0,0];  
    case 'BM'
          border = 0;
         %offset = [0,-151];
         offset = [0,0];
         
    case 'MBM'
       border = 0;
     offset = [0,0];  
end

fvert = gi_flat.vertices(:,1:2);
lower = min(fvert,[],1);
lower(:) = min(lower);
upper = max(fvert,[],1);
upper(:) = max(upper);
fvert = fvert - repmat(lower,size(fvert,1),1);
fvert = fvert ./ repmat(upper-lower,size(fvert,1),1);
fvert(:,1) = fvert(:,1)*(f_shape(1)-border-1)+border/2+offset(1);
fvert(:,2) = fvert(:,2)*(f_shape(2)-border-1)+border/2+offset(2);
fm = zeros(f_shape);

switch FM
    case 'BM'
         indx = sub2ind(f_shape,f_shape(1)-round(fvert(:,1)),round(fvert(:,2)));
    case 'MBM'
        T = eye(3);
        T(1:2,3) = -f_shape/2;
        T2 = T;
        T2(1:2,3) = f_shape/2-[-50,90]; 
        
        S = eye(3);
        S(1,1) = 1.2;
        S(2,2) = 1.2;

        R = eye(3);
        %alpha = 30 * 2*pi/360;
        
        alpha = 17 * 2*pi/360;
        S(1,1) = 1.15;
        S(2,2) = 1.15;
        T2(1:2,3) = f_shape/2-[-70,70]; 
        
        R(1:2,1:2) = [cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
        M = T2*R*S*T;
        verts_ = cat(2,fvert,ones([size(fvert,1),1]));
        fvert_ = (M*verts_')'; 
         indx = sub2ind(f_shape,f_shape(2)-round(fvert_(:,2)),f_shape(1)-round(fvert_(:,1)));
         fvert = verts_(:,1:3);
end
fm(indx(valid)) = 1;
fm(indx(~valid)) = -1;
valid = fm(indx)>0;

figure(1);imagesc(fm)
valid_sample_pts = valid;
imwrite(cat(3,(fm<0)*255,(fm>0)*255,(fm<0)*0),[ofolder,'/flatmap_samples.png']);
mat_2D = fvert(:,1:2);
mat_indx = indx;
save([ofolder,'/verts_in_stack_location.mat'],'mat_2D','mat_indx','valid_sample_pts')
return
%%
fm2 = zeros(size(fm));
fm2(indx(valid))=1;
assert(sum((fm2(:)+(fm(:)<0))==2)==0);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  flatmap valid mask computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FM
    case 'Rosa'
        fm_remove_mask = imclose(fm==-1,strel('disk',5));
        fm_mask = imclose((fm==1)&(~fm_remove_mask),strel('disk',15));
        fm_mask = imopen(fm_mask,strel('disk',5));
        figure(1);imagesc(cat(3,fm_mask,fm*0.5+0.5,fm*0.5+0.5))        
    case 'BM'
        
        fm_remove_mask = imclose(fm==-1,strel('disk',5));
        fm_mask = imclose((fm==1)&(~fm_remove_mask),strel('disk',35));
        fm_mask = imopen(fm_mask,strel('disk',5));
        
        if exist('fm_mask_2','var') 
            fm_mask = fm_mask & fm_mask_2;
        end
        figure(1);imagesc(cat(3,fm_mask,fm*0.5+0.5,fm*0.5+0.5))
    case 'MBM'

        if highres
        fm_remove_mask = imdilate(fm==-1,strel('disk',3));
        fm_mask = imclose((fm==1)&(~fm_remove_mask),strel('disk',15));
        fm_mask = imopen(fm_mask,strel('disk',5));
        figure(1);imagesc(cat(3,fm_mask,fm*0.5+0.5,fm*0.5+0.5))           
        else
        fm_remove_mask = imclose(fm==-1,strel('disk',30));
        fm_mask = imclose((fm==1)&(~fm_remove_mask),strel('disk',6));
        fm_mask = imopen(fm_mask,strel('disk',10));
        figure(1);imagesc(cat(3,fm_mask,fm*0.5+0.5,fm*0.5+0.5))   
        end
end

imwrite(cat(3,fm_mask,fm*0.5+0.5,fm*0.5+0.5),[ofolder,'/flatmap_mask.png'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  fm interpolation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = ndgrid(1:f_shape(1),1:f_shape(2));
P = cat(2,double(X(indx(valid))),double(Y(indx(valid))));
V = double(indx(valid));
[V_unique,P_unique] = groupsummary(V,P,@max);
P_unique = [P_unique{:}];
I = scatteredInterpolant(P_unique,V_unique,'nearest');
I.ExtrapolationMethod = 'none';

fm_indx_interp = zeros(f_shape);
fm_indx_interp(:)=I(X(:),Y(:));
fm_indx_interp(isnan(fm_indx_interp(:))) = 0;
fm_indx_interp_all = fm_indx_interp ;
fm_indx_interp = fm_indx_interp .* fm_mask;

figure(1);imagesc(fm_indx_interp);

%%
if ~highres
     [X,Y] = ind2sub(size(fm_indx_interp_all),indx(:));
     save([ofolder,'/flatmap2flatmapstack',FM,'.mat'],'X','Y');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  create flatmapstack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%
f_shape = [1000,1000];
fm_stack_indx = zeros([f_shape,smpl_pts,3],'uint32');
for d = 1:smpl_pts
    fprintf('%d %d\n',d,smpl_pts);
    pts = squeeze(result_tracks_array(:,d,:))+1;
    for a = 1:3
        tmp = zeros(f_shape,'uint32');
        pts_ = round(pts(:,a));
        tmp(indx(valid_sample_pts(:))) = pts_(valid_sample_pts(:));
        tmp(fm_indx_interp(:)>0) = tmp(fm_indx_interp(fm_indx_interp(:)>0));
        tmp(~fm_mask) = 0;
        fm_stack_indx(:,:,d,a) = tmp;
    end
end

fm_stack_indx_xyz = fm_stack_indx;
%%
nifti_stack = make_nii(fm_stack_indx);
save_nii(nifti_stack,[ofolder,'/fm_stack_XYZ_1000_TC_std.nii.gz']);



%%
f_shape = [1000,1000];
fm_stack_indx = zeros([f_shape,smpl_pts],'uint32');
for d = 1:smpl_pts
    fprintf('%d %d\n',d,smpl_pts);
    pts = squeeze(result_tracks_array(:,d,:))+1;
    indx3D = sub2ind(shape,round(pts(:,1)),round(pts(:,2)),round(pts(:,3)));
    tmp = zeros(f_shape,'uint32');
    
    tmp(indx(valid(:))) = indx3D(valid(:));
    tmp(fm_indx_interp(:)>0) = tmp(fm_indx_interp(fm_indx_interp(:)>0));
    tmp(~fm_mask) = 0;
    fm_stack_indx(:,:,d) = tmp;
end

   
%%
nifti_stack = make_nii(fm_stack_indx);
save_nii(nifti_stack,[ofolder,'/fm_stack_index_1000_TC_std.nii.gz']);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DEFORMATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~highres

   
    valid = max(squeeze(result_tracks_array(:,:,1)),[],2) > 0;
    f_stack_shape = [500,500];
    fm_stack_deform = zeros([f_stack_shape,smpl_pts],'single');
    fm_stack_dist = zeros([f_stack_shape,smpl_pts],'single');
    for d = 1:smpl_pts
        fprintf('%d %d\n',d,smpl_pts);
        tmp = zeros(f_shape,'single');

        vertex_feat_ = squeeze(vertex_feat(:,d));
        tmp(indx(valid(:))) = (vertex_feat_(valid(:)));
        tmp(fm_indx_interp(:)>0) = tmp(fm_indx_interp(fm_indx_interp(:)>0));

       
        fm_stack_deform(:,:,d) = imresize(tmp,f_stack_shape,'bilinear');
        
        tmp = zeros(f_shape,'single');
        tmp(indx(valid(:))) = DIc(:,d);
        tmp(fm_indx_interp(:)>0) = tmp(fm_indx_interp(fm_indx_interp(:)>0));
        fm_stack_dist(:,:,d) = imresize(tmp,f_stack_shape,'bilinear');
    end

    %%
    nifti_stack = make_nii(fm_stack_deform);
    save_nii(nifti_stack,[ofolder,'/fm_stack_deform.nii.gz']);
    %%
    
    nifti_stack = make_nii(fm_stack_dist);
    save_nii(nifti_stack,[ofolder,'/fm_stack_dist.nii.gz']);
    
    %%
    TC_img_deformation = zeros(shape);
    fm_stack_deform_ =  zeros([f_shape,smpl_pts],'single');
    for d = 1:smpl_pts
        fm_stack_deform_(:,:,d) = imresize(squeeze(fm_stack_deform(:,:,d)),[1000,1000]);
    end
    TC_img_deformation(fm_stack_indx(fm_stack_indx(:)>0)) = fm_stack_deform_(fm_stack_indx(:)>0);
    %%
    [X,Y,Z] = ndgrid(1:shape(1),1:shape(2),1:shape(3));
    valid_def = TC_img_deformation>0;
    P = cat(2,double(X(valid_def)),double(Y(valid_def)),double(Z(valid_def)));


    I = scatteredInterpolant(P,TC_img_deformation(valid_def(:)),'linear');
    I.ExtrapolationMethod = 'none';
    %%
    valid_def = label2>0;
    P_ = cat(2,double(X(valid_def)),double(Y(valid_def)),double(Z(valid_def)));
    TC_img_deformation(valid_def(:)) = I(P_);
    
    %%

    trace_nii = TC_img;
    trace_nii.hdr.dime.datatype = 16;
    trace_nii.hdr.dime.bitpix = 32;
    trace_nii.hdr.dime.glmax = max(TC_img_deformation(:));
    trace_nii.hdr.dime.glmin = min(TC_img_deformation(:));
    trace_nii.img = single(TC_img_deformation);
    save_untouch_nii(trace_nii,[ofolder,'/deformations_TC_std.nii.gz']);

end


