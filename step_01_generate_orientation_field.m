
%%
label2 = load_untouch_nii('data/cortex_boundary_TC_std.nii.gz');
label2 = label2.img;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   order 0 orientation field smoothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mask = label2;
mask = (mask == 2) | (mask == 3);
data =  zeros(size(label2));
data(:)=0.5;
data(label2(:) == 2) = 0;
data(label2(:) == 3) = 1;

H = zeros([3,3,3],'double');
    H(2,2,1) = 1;
    H(:,:,2) = [0,1,0;1,-6,1;0,1,0];
    H(2,2,3) = 1;
%%
iters = 5000;


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

cortex_dist = single(gather(data2));



%%

trace_nii = label2;
trace_nii.hdr.dime.datatype = 16;
trace_nii.hdr.dime.bitpix = 32;
trace_nii.img = single(cortex_dist);
save_untouch_nii(trace_nii,'cortex_dist.nii.gz');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   order 1 orientation field smoothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
gm = stafield(single(label2 == 1));
kernel = stafield('gauss',size(label2),1.5,0,0,'STA_FIELD_STORAGE_R','single');
gms = gm.fft().prod(kernel.fft(),0).ifft();
gmd = gms.deriv(1);

%%
gmd_n = gmd.norm();
gmd_n = 1/(squeeze(gmd_n.data(1,1,:,:,:))+eps);
gmdn = gmd;
for a = 1:2
    for r = 1:2
        gmdn.data(a,r,:,:,:) = squeeze(gmdn.data(a,r,:,:,:))  .* gmd_n;
    end
end
%start_VIannotator(squeeze(gmdn.data(1,1,:,:,:)));
%%
step = 1;
mask = label2;
mask = (mask == 2) | (mask == 3);
gm_in = gmdn;
gm_out = gmdn;
for m=1:2
    for r = 1:2
        gm_in.data(m,r,:,:,:) = squeeze(gmdn.data(m,r,:,:,:)).*(label2 == 2);
        gm_out.data(m,r,:,:,:) = -squeeze(gmdn.data(m,r,:,:,:)).*(label2 == 3);
    end
end
gm_both = gm_in;
gm_both.data = gm_in.data + gm_out.data; 
%%


result = gm_both;
%result.data = 
%%
iters = 5000;
for m = 1:size(gm_both.data,2)
    for r=1:2
        fprintf(' %d %d  \n',m,r);
        data = squeeze(gm_both.data(r,m,1:step:end,1:step:end,1:step:end));

        data2_org = gpuArray(data);

        data2 = gpuArray(data);
        mask2 = gpuArray(mask);
        H2 =  gpuArray(H);
        if sum(abs(data(:)))>0.000001
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
        else
            fprintf('skipping %d %d  \n',m,r);
        end
        result.data(r,m,:,:,:) = gather(data2);
      
    end
end
%%
save('data/result_t1_n.mat','result','-v7.3')
gradf = sta_s2c(result);
tracer_mask = single(label2 >0);
tracer_gradf = gradf;
save('data/fmstack_data_n.mat','tracer_mask','tracer_gradf','cortex_dist','-v7.3')

