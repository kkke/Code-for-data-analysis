clear
gcp;
%% registration
name = 'F:\Imaging in GC\RVKC340\061218\RVKC340_180612_5avg.tif';
tic; Y = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)
%%
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));
%Y = Y - min(Y(:));
% set parameters (first try out rigid motion correction)
%% rigid motion correction 
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
% perform motion correction
tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
%%
%% now try non-rigid motion correction (also in parallel)
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc

%% for displaying the movie
% ref = Y(:,:,1);
% nnY = quantile(ref(:),0.005);
% mmY = quantile(ref(:),0.995);
% 
% %% plot a movie with the results
% figure;
% for t = 1:1:T
%     subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(M1(:,:,t),[nnY,mmY]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     drawnow;
%     pause(0.02);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write the tiff file
% M1 = int16(M1);
res = saveastiff(M1, 'reg1.tif')
res = saveastiff(M2, 'reg2.tif')