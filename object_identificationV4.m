%% Select Image
[fileX,pathX] = uigetfile('*.tif');
disp("Opening File:" + fileX);
%% Load Image 
%8-15-2018
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_8_15 MSM HEK DC-EF/8_15_2018 HEK DC-EF MSM - No stim_2/8_15_2018 HEK DC-EF MSM - No stim_2_MMStack_Pos0.ome.tif';
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_8_15 MSM HEK DC-EF/8_15_2018 HEK DC-EF MSM - stim 20V every 10s_1/8_15_2018 HEK DC-EF MSM - stim 20V every 10s_1_MMStack_Pos0.ome.tif';
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_8_15 MSM HEK DC-EF/8_15_2018 HEK DC-EF MSM - stim 10V every 10s_1/8_15_2018 HEK DC-EF MSM - stim 10V every 10s_1_MMStack_Pos0.ome.tif';

%8-17-2018
    %IN MEDIA
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in media_1/2018_8_17 - HEK-NK _ in media_1_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in media_2/2018_8_17 - HEK-NK _ in media_2_MMStack_Pos0.ome.tif';
    %IN ACSF
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_2/2018_8_17 - HEK-NK _ in ACSF_2_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_3/2018_8_17 - HEK-NK _ in ACSF_3_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_4/2018_8_17 - HEK-NK _ in ACSF_4_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_5/2018_8_17 - HEK-NK _ in ACSF_5_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_6/2018_8_17 - HEK-NK _ in ACSF_6_MMStack_Pos0.ome.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_17 HEK-NK in Media followed by ACSF/2018_8_17 - HEK-NK _ in ACSF_7/2018_8_17 - HEK-NK _ in ACSF_7_MMStack_Pos0.ome.tif';
       
%Older Vids of HEK
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Losert Matlab Testing/Picking up Cytoplasm features/Test_Images/2018_5_21_HEK-NK_DC-EF_with_CellLight-488BeRST640.mvd2 - 60x_stimulation every 30s_10V field_640 and 488 _near bottom_series 21.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Losert Matlab Testing/Picking up Cytoplasm features/Test_Images/2018_7_16_HEK DC EF P14 in media only_P33 and 14 in ACSF.mvd2 - 40x_P14_in ACSF_12.5V stim at 1,130,2,230_3.5min 3.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Losert Matlab Testing/Picking up Cytoplasm features/Test_Images/2018_6_15_HEK-NK_Cellight+BeRST_Lowest voltage_and_hysteresis.mvd2 - 60x_10V every 30s for 6min_then15V every 30s for 6min_12mintotal.tif';
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Losert Matlab Testing/Picking up Cytoplasm features/Test_Images/2018_5_21_HEK-NK_DC-EF_with_CellLight-488BeRST640.mvd2 - 60x_stimulation every 30s_10V field_640 and 488 - 3min 1.5 timepointss-1_BeRSTonly.tif';
        
%Neurons
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Spinning Disk/2018-07-27 Hc neurons from 2018-07-20 + CellLight GFP + BeRST [From KMO]/ImageJ extractions/2018-07-27 Hc neurons from 2018-07-20 + CellLight GFP + BeRST.mvd2 - B2_40x_timelapse3_Actin1x+BeRST8x_5min (All 8 frames sequentially).tif';
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Spinning Disk/2018-08-01 Cx neurons from 2018-07-30 + BeRST [From KMO]/imagej extractions/2018-08-01 Cx neurons from 2018-07-30 + BeRST.mvd2 - dish1_BeRST_1-5-1_DC_movie2_minutue 1-2.tif';    
    
    
%8-15-2018
    %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_15 MSM HEK DC-EF/8_15_2018 HEK DC-EF MSM - No stim_2/8_15_2018 HEK DC-EF MSM - No stim_2_MMStack_Pos0.ome.tif';
    
%8-20-2018 (all have 14ms interval ontop of the limit of 13.8ms exposure)
    %In ACSF
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_20 HEK-NK in Media followed by ACSF/2018_8_20 - HEK-NK _ in ACSF_8/2018_8_20 - HEK-NK _ in ACSF_8_MMStack_Pos0.ome.tif';
    %in Media
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_20 HEK-NK in Media followed by ACSF/2018_8_20 - HEK-NK _ in Media_9/2018_8_20 - HEK-NK _ in Media_9_MMStack_Pos0.ome.tif';
        
%8-22-2018
    %in Plate 2 ACSF 
        %im_path = '/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_22/2018_8_22 - HEK-NK _Plate2_ACSF(with 5ms interval)_1/2018_8_20 - HEK-NK _Plate2_ACSF(with 5ms interval)_2_MMStack_Pos0.ome.tif';
        im_path = fullfile(pathX,fileX);
%% Preloaded Information 

% Information from Tiff 
info = imfinfo(im_path);

    %Height and Width of original video 
    im_height = info(1).Height;
    im_width = info(1).Width;
    
    
    %Frame Range of original video
    max_no_frame = numel(info);
    max_f = length(info);
    First_frame = 1;           % <--- First frame you want to quantify
    Last_frame = max_no_frame;   
    
N = 0;
Y = 1;
Ones = ones(im_height,im_width);

%Filters 
    filt4 =[-1 -1 -10 -1 -1 ; -1 -10 -100 -10 -1 ; -10 -100 491 -100 -10 ; -1 -10 -100 -10 -1 ; -1 -1 -10 -1 -1];

    filtsize = 3;
    filt3 = ones(filtsize) .* (-1) ;
    midindx = (round(filtsize/2));
    filt3(midindx ,midindx) = (filtsize*filtsize-1);  
%% User input - Preset Analysis Parameters 

%fs = 13.89+14;  %the freq at which the imgs are taken (eg: 1000/[frame per ms]
%fs = round(1000/fs); %sampling frequency, in Hz(eg: frames per second);
%fs = max_no_frame/(60*6)
fs = 15.3;

%Analysis Parameters
F_RANGE = 1:max_no_frame;
C_Analysis_Percent = 100; %In percentage, the number of cells analyzed 


USE_TRACKING = 0;  % value of 0 equals off, % value of 1 equals on

%Buffer for ROIs
dilate_size = 20;

%IF looking at small cells, use skeleton and set to 1
USE_Skel = 0;

%Modifiers to identify boundaries
Thresh_Mod = 0.3;       %Adjust until BG and cells are separated cleanly 
Clean_Mod = 1.2;
Purne_Extensions_Mod = 50;

%Filter Size for objects
Min_Obj_Size = 50;
Max_Obj_Size = 800;


%Output Paramters
line_style = '-';
axis_y_lim = 0.4;

indvplot_x = 5;
indvplot_y = 5;

Resize_factor = 1;

%% Load Original Into Storage 
% 
% %start stopclock
% stopclock = 0;
% tic;
% 
% Storage = struct([]);
% %prog_bar = waitbar(0,'Loading Frames Into Storage');
% for frame = F_RANGE
%     Storage(frame).OG_im = imread(im_path,frame); 
%      %frame_percent = frame/length(F_RANGE);
%      %waitbar(frame_percent ,prog_bar, sprintf("Loading Frames "+F_RANGE(1)+ "-"+max(F_RANGE)+ "\n" +round(frame_percent.*100)+"%% Completed"));
% end
% %close(prog_bar);
% toc;
% stopclock = stopclock + toc;
%% Load Original Into Storage (Much Faster)
%start stopclock
stopclock = 0;
tic;
warning ('off','all');
FileTif= im_path; %fullfile(pathX,fileX);
%InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
%FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
Storage = struct([]);
 
TifLink = Tiff(FileTif, 'r');
for frame=F_RANGE
   TifLink.setDirectory(frame);
   Storage(frame).OG_im = imresize((TifLink.read()),Resize_factor);
   %FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();
toc;
stopclock = stopclock + toc;
%% Segmentation and Storage/Overview 
tic;


if USE_TRACKING == 1
    F_RANGE_ROI = F_RANGE;
else
    F_RANGE_ROI = F_RANGE(1);
end

cont = 0;
while cont == 0
SEclose = strel('disk',2);  % vlaue of 2 seems to work best
SE_bg = strel('disk',round (sqrt(im_width * im_height)*0.1));
SEerode = strel('disk',dilate_size);
SEdilate = strel('disk',dilate_size);


%Create ROI/Cell Masks
for frame = F_RANGE_ROI
    
    %load in image
    im = Storage(frame).OG_im;
    im = im2double(im);

    %adjust image and even out background gradient
    im0 = imadjust(im);
    im0 = imgaussfilt(im0,Clean_Mod);
    im0 = im0 - imopen(im0,SE_bg); %even out bg gradient so that it can be picked up 
    level = graythresh(im0); % calculate level using Otsu method to try and separate 
    mask = im0>level .* Thresh_Mod; % mask that tries to select only Objects and remove dark BG

    %BG_ROI = bwareafilt(logical(mask<1),[Max_Obj_Size*1.1 (10^120)]);
    BG_ROI = logical(mask<1);
    BG_ROI = imerode(BG_ROI,SEerode);
    
    %Filter image to better see border
    im_Filt = im0;
    im_Filt = filter2 (filt3,im0); % runs image through filter to find borders
    im_Filt = imadjust(imgaussfilt(double(im_Filt),Clean_Mod)); %blurs to remove noise [val of 1-2 works best]
    im_Filt = imclose(im_Filt,SEclose);
    
    %skeltonize borders to segment image
    im_Skel = logical(bwmorph(im_Filt,'thin',inf) .* mask);
    im_Skel = bwareafilt(im_Skel,[Min_Obj_Size 10^120000]);
    im_Skel = bwmorph(im_Skel,'spur',Purne_Extensions_Mod);
    im_Skel = imgaussfilt(double(im_Skel),1);
    
    %invert image
    im_Inv = im_Skel==0;

    %remove particles outiside of range 
    im_Cull = bwareafilt(logical(im_Inv), [Min_Obj_Size Max_Obj_Size]);
    
    %create labeled image and acquire number of unqiue cells/obj picked up
    labl = bwlabel(im_Cull);
    num_labl = max(labl(:));

    Storage(frame).Objs = labl;
    Storage(frame).Objs_Num = num_labl;
    
%Preview and Checkpoint at 1st frame
    if frame == F_RANGE_ROI
         fig(1) = figure('Name','Overview of ROI Creation'); 
         subplot(331); imagesc(im); axis 'image'; title("Initial Image");  
         subplot(332); imagesc(im0); axis 'image'; title("ImAdjusted"); 
         subplot(333); imagesc(mask); axis 'image'; title("Separating ROI from BG"); 
         subplot(334); imagesc(im_Filt); axis 'image'; title("Filtered and ImAdjusted"); 
         subplot(335); imagesc(im_Skel); axis 'image'; title("Segmentation");
         subplot(336); imagesc(im_Inv); axis 'image'; title("Inverted"); 
         subplot(337); imagesc(im_Cull); axis 'image'; title("Culled large Areas");
         subplot(338); imagesc(labl); axis 'image'; title(["Cells/Objects Identified: "+num_labl," Over "+length(F_RANGE)+" frames"]);
         subplot(339); imagesc(imdilate(im_Cull, SEdilate) - im_Cull); axis 'image'; title("Dilated Edges"); 
 
         fuseROI = imfuse(imadjust(im),imdilate(im_Cull, SEdilate) - im_Cull);
         fuseBG_ROI = imfuse(imadjust(im),BG_ROI);
 
         fig(2) = figure('Name','Overview of ROIs'); 
         subplot(121); imagesc(fuseROI); axis 'image'; title("Cell/Obj ROI Preview"); 
         subplot(122); imagesc(fuseBG_ROI); axis 'image'; title("BG ROI Preview");

         movegui(fig(1), 'northeast');
         movegui(fig(2), 'northwest');
         
        prompt1 = char(["Proceed to Analysis?";"0 = No"; "1 = Yes"]);
        title1 = "Checkpoint";
        dims1 = [1 35]; 
        checkpoint_dlg = inputdlg(prompt1,title1,dims1);
        continue_var = str2double(checkpoint_dlg);

    
        %Checkpoint allows to continue or stop and change values
         if continue_var == 0
             prompt2 = {'Cell Analysis Percentage','Sampling Frequency in Hz (frames per second)','Modify ROI Dilation Size',sprintf('Modify BG Thresh [0-1] \nHigher values increase coverage \nLower values decrease coverage'),sprintf('Modify Cleaning Value [>1] \nHigher values clean image more \nLower values sharpen image more'),'Modify Pruning Value','Min Cell/Obj Area(px)','Max Cell/Obj Area(px)'};
             title2 = 'Input';
             dims2 = [1 35];
             definput2 = {num2str(C_Analysis_Percent),num2str(fs),num2str(dilate_size),num2str(Thresh_Mod),num2str(Clean_Mod),num2str(Purne_Extensions_Mod),num2str(Min_Obj_Size),num2str(Max_Obj_Size)};
             answers_vals = inputdlg(prompt2,title2,dims2,definput2);
             
             C_Analysis_Percent = str2double(answers_vals(1));
             fs = str2double(answers_vals(2));
             dilate_size = str2double(answers_vals(3));
             Thresh_Mod = str2double(answers_vals(4));
             Clean_Mod = str2double(answers_vals(5));
             Purne_Extensions_Mod =  str2double(answers_vals(6));
             Min_Obj_Size = str2double(answers_vals(7));
             Max_Obj_Size = str2double(answers_vals(8));
             
             close all;
% %            return
         elseif continue_var == 1
             disp("Continuing Onto Processing Frames");
             cont = 1;
             close all;
         end    
    pause(0.00001);
    end
    
end

end
%% Segmentation and Storage/Overview For Neurons 
Use_Neuron = 0;

if Use_Neuron ==1; 
tic;

if USE_TRACKING == 1
    F_RANGE_ROI = F_RANGE;
else
    F_RANGE_ROI = F_RANGE(1);
end

cont = 0;
while cont == 0
    
SEclose = strel('disk',2);  % vlaue of 2 seems to work best
SE_bg = strel('disk',round (sqrt(im_width * im_height)*0.1));
SEerode = strel('disk',dilate_size);
SEdilate = strel('disk',dilate_size);


%Create ROI/Cell Masks
for frame = F_RANGE_ROI
    
    im = Storage(frame).OG_im;
    im = im2double(im);
    
    im0 = im; 
    filt = bwareafilt(im0<0.06,[20 110]);
    
    labl = bwlabel(filt);
    num_labl = max(labl(:));

    Storage(frame).Objs = labl;
    Storage(frame).Objs_Num = num_labl;
    cont = cont+1;
%Preview and Checkpoint at 1st frame
%     if frame == F_RANGE(1)
%          fig(1) = figure('Name','Overview of ROI Creation'); 
%          subplot(331); imagesc(im); axis 'image'; title("Initial Image");  
%          subplot(332); imagesc(im0); axis 'image'; title("ImAdjusted"); 
%          subplot(333); imagesc(mask); axis 'image'; title("Separating ROI from BG"); 
%          subplot(334); imagesc(im_Filt); axis 'image'; title("Filtered and ImAdjusted"); 
%          subplot(335); imagesc(im_Skel); axis 'image'; title("Segmentation");
%          subplot(336); imagesc(im_Inv); axis 'image'; title("Inverted"); 
%          subplot(337); imagesc(im_Cull); axis 'image'; title("Culled large Areas");
%          subplot(338); imagesc(labl); axis 'image'; title(["Cells/Objects Identified: "+num_labl," Over "+length(F_RANGE)+" frames"]);
%          subplot(339); imagesc(imdilate(im_Cull, SEdilate) - im_Cull); axis 'image'; title("Dilated Edges"); 
%  
%          fuseROI = imfuse(imadjust(im),imdilate(im_Cull, SEdilate) - im_Cull);
%          fuseBG_ROI = imfuse(imadjust(im),BG_ROI);
%  
%          fig(2) = figure('Name','Overview of ROIs'); 
%          subplot(121); imagesc(fuseROI); axis 'image'; title("Cell/Obj ROI Preview"); 
%          subplot(122); imagesc(fuseBG_ROI); axis 'image'; title("BG ROI Preview");
% 
%          movegui(fig(1), 'northeast');
%          movegui(fig(2), 'northwest');
%          
%         prompt1 = char(["Proceed to Analysis?";"0 = No"; "1 = Yes"]);
%         title1 = "Checkpoint";
%         dims1 = [1 35]; 
%         checkpoint_dlg = inputdlg(prompt1,title1,dims1);
%         continue_var = str2double(checkpoint_dlg);
% 
%     
%         %Checkpoint allows to continue or stop and change values
%          if continue_var == 0
%              prompt2 = {'Cell Analysis Percentage','Sampling Frequency in Hz (frames per second)','Modify ROI Dilation Size',sprintf('Modify BG Thresh [0-1] \nHigher values increase coverage \nLower values decrease coverage'),sprintf('Modify Cleaning Value [>1] \nHigher values clean image more \nLower values sharpen image more'),'Modify Pruning Value','Min Cell/Obj Area(px)','Max Cell/Obj Area(px)'};
%              title2 = 'Input';
%              dims2 = [1 35];
%              definput2 = {num2str(C_Analysis_Percent),num2str(fs),num2str(dilate_size),num2str(Thresh_Mod),num2str(Clean_Mod),num2str(Purne_Extensions_Mod),num2str(Min_Obj_Size),num2str(Max_Obj_Size)};
%              answers_vals = inputdlg(prompt2,title2,dims2,definput2);
%              
%              C_Analysis_Percent = str2double(answers_vals(1));
%              fs = str2double(answers_vals(2));
%              dilate_size = str2double(answers_vals(3));
%              Thresh_Mod = str2double(answers_vals(4));
%              Clean_Mod = str2double(answers_vals(5));
%              Purne_Extensions_Mod =  str2double(answers_vals(6));
%              Min_Obj_Size = str2double(answers_vals(7));
%              Max_Obj_Size = str2double(answers_vals(8));
%              
%              close all;
% % %            return
%          elseif continue_var == 1
%              disp("Continuing Onto Processing Frames");
%              cont = 1;
%              close all;
%          end    
%     pause(0.00001);
%     end
    
end
 Max_BG_ROI = imerode(filt,SEerode);
end


end
%% BG ROI Creation
%create MAX BG ROI
im_max = zeros(im_height*Resize_factor,im_width*Resize_factor);
Max_BG_ROI = ones(im_height*Resize_factor,im_width*Resize_factor);

for frame = F_RANGE
    im = Storage(frame).OG_im;
    im = im2double(im);

    %adjust image and even out background gradient
    im0 = imadjust(im);
    im0 = imgaussfilt(im0,Clean_Mod);
    im_max = max(im_max, im0); 
    
end

    im1 = im_max - imopen(im_max,SE_bg); %even out bg gradient so that it can be picked up 
    level = graythresh(im1); % calculate level using Otsu method to try and separate 
    mask = im1>level .* Thresh_Mod; % mask that tries to select only Objects and remove dark BG

    %BG_ROI = bwareafilt(logical(mask<1),[Max_Obj_Size*1.1 (10^120)]);
    BG_ROI = logical(mask<1);
    BG_ROI = imerode(BG_ROI,SEerode);
    Max_BG_ROI = Max_BG_ROI .* BG_ROI;

toc;
stopclock = stopclock + toc;
%% Final Preview 
close all;
fig(1) = figure('Name','Overview of ROI Creation');
subplot(331); imagesc(im); axis 'image'; title("Initial Image");
subplot(332); imagesc(im0); axis 'image'; title("ImAdjusted"); 
subplot(333); imagesc(mask); axis 'image'; title("Separating ROI from BG");
subplot(334); imagesc(im_Filt); axis 'image'; title("Filtered and ImAdjusted"); 
subplot(335); imagesc(im_Skel); axis 'image'; title("Segmentation");
subplot(336); imagesc(im_Inv); axis 'image'; title("Inverted"); 
subplot(337); imagesc(im_Cull); axis 'image'; title("Culled large Areas"); 
subplot(338); imagesc(labl); axis 'image'; title("Cells/Objects Identified: "+num_labl);
subplot(339); imagesc(imdilate(im_Cull, SEdilate) - im_Cull); axis 'image'; title("Dilated Edges"); 

fuseROI = imfuse(imadjust(im),imdilate(im_Cull, SEdilate) - im_Cull);
fuseBG_ROI = imfuse(imadjust(im),Max_BG_ROI);

fig(2) = figure('Name','Overview of ROIs'); 
subplot(121); imagesc(fuseROI); axis 'image'; title("Cell/Obj ROI Preview"); 
subplot(122); imagesc(fuseBG_ROI); axis 'image'; title("Max BG ROI Preview"); 

movegui(fig(1), 'northeast');
movegui(fig(2), 'northwest');

pause(0.00001);
%% Cell_Analysis Range
%this will do a certain percentage of total cells found
C_RANGE = 1:round(Storage(F_RANGE(1)).Objs_Num/(Storage(F_RANGE(1)).Objs_Num * (C_Analysis_Percent/100) )):Storage(F_RANGE(1)).Objs_Num;

%use this if using or analyzing specific cells
%C_RANGE = 1:Storage(F_RANGE(1)).Objs_Num;
%% Gather BG Raw Data

Data = struct([]);
Data(1).BG_Raw_IntArea = [];
%gather BG data
tic;
for frame =  F_RANGE
    mult = Max_BG_ROI .* double(Storage(frame).OG_im);
    Data(1).BG_Raw_IntArea(frame,1) = sum(mult(:)) / sum(Max_BG_ROI(:));
end
toc;
stopclock = stopclock + toc;
%% Gather Membrane Raw Data
%gather cell data
tic;



% calculte the cell mask first then run through each frame 
prog_bar = waitbar(0,'Gathering Cell/Obj ROI For Each Frame');
Data.Cell_Raw_IntArea = [];
for frame = F_RANGE
    for cell = C_RANGE
        roi = (imdilate(Storage(F_RANGE(1)).Objs==cell, SEdilate) - (Storage(F_RANGE(1)).Objs==cell));
        mult = roi .* double(Storage(frame).OG_im); 
        Data(1).Cell_Raw_IntArea(frame,cell) = sum(mult(:)) / sum(roi(:)) ;
    end
    frame_percent = (frame-F_RANGE(1)+1)/length(F_RANGE);
    waitbar(frame_percent ,prog_bar, sprintf("Gathering Cell/Obj ROI\n"+round(frame_percent.*100)+"%% Completed"));
end
close(prog_bar);
toc;
stopclock = stopclock + toc;
%% Gather Cell Metrics

for frame = F_RANGE_ROI
    for cell = C_RANGE
    stats(cell) = regionprops(Storage(frame).Objs==cell,'Area','Perimeter','Eccentricity','Orientation');
    end
end
%% Calculate True Int and DF/F Based on Data
window_size = ceil(length(F_RANGE)*0.20);
w_sz = floor(window_size/2); 

%calculate True int
tic;
Data.Cell_TrueInt = [];
for frame = F_RANGE
    Data(1).Cell_TrueInt(frame,:) = Data(1).Cell_Raw_IntArea(frame,:) - Data(1).BG_Raw_IntArea(frame,1);
end

%calculate Delta F of F for BG
Data.BG_DFoF = [];
for frame = F_RANGE
    
    leftbound = frame - w_sz; 
    rightbound = frame + w_sz;
    windowRng = leftbound:rightbound;
    if sum(max(windowRng)==F_RANGE)==0 && sum(min(windowRng)==F_RANGE)==0
        windowRng = windowRng;
    elseif sum(max(windowRng)==F_RANGE)==1 && sum(min(windowRng)==F_RANGE)==0
        windowRng =  F_RANGE(1) : frame + ((2*w_sz) - (frame-F_RANGE(1)));
    elseif sum(max(windowRng)==F_RANGE)==0 && sum(min(windowRng)==F_RANGE)==1
        windowRng = frame - ((2*w_sz)-(max(F_RANGE) - frame)) : max(F_RANGE);
    end
    
    currentInt = Data.BG_Raw_IntArea(frame);
    initialInt = Data.BG_Raw_IntArea(F_RANGE(1));
%    initialInt = sum(Data(1).Cell_TrueInt(window,cell)) / length(window);
    Data(1).BG_DFoF(frame,1) = (currentInt-initialInt)/(initialInt);
end


%calculate Delta F of F for Cells
Data.Cell_DFoF = [];
for frame = F_RANGE
    
    %Use if using Window
    leftbound = frame - w_sz; 
    rightbound = frame + w_sz;
    windowRng =  leftbound:rightbound;
    if sum(max(windowRng)==F_RANGE)==0 && sum(min(windowRng)==F_RANGE)==0
        %IF the max and min of window ranges are within the Frame range then leave window range the same 
        windowRng = windowRng;
    elseif sum(max(windowRng)==F_RANGE)==1 && sum(min(windowRng)==F_RANGE)==0
        %If the min is outside the window range
        windowRng =  F_RANGE(1) : frame + ((2*w_sz) - (frame-F_RANGE(1))) ;
    elseif sum(max(windowRng)==F_RANGE)==0 && sum(min(windowRng)==F_RANGE)==1
        windowRng = frame - ((2*w_sz)-(max(F_RANGE) - frame)) : max(F_RANGE);
    end
    
    for cell = C_RANGE
        currentInt = Data(1).Cell_TrueInt(frame,cell);
        
        %based on first frame
        %initialInt = Data(1).Cell_TrueInt(F_RANGE(1),cell);   
        
        %based on mean/average
        %initialInt = mean(Data(1).Cell_TrueInt(window,cell));   
        
        %based on minimum
        %initialInt = mean(mink(Data(1).Cell_TrueInt(window,cell),2)) ; 
        
%         middle10_start = round(size(Data(1).Cell_TrueInt,1)*.5) - round((size(Data(1).Cell_TrueInt,1)*.01)/2);
%         middle10_end = middle10_start + round((size(Data(1).Cell_TrueInt,1)*.1));
%      
%         sorted_cell = sort(Data(1).Cell_TrueInt(:,cell));
%         sorted_cell_mid = sorted_cell(middle10_start:middle10_end);
%         initialInt = mean(Data.sorted_cell_mid); 
        
        %based on percentile in the [] - usually 25, but 10 if decay is severe
        initialInt = prctile(Data(1).Cell_TrueInt(windowRng,cell),10,1); 

        

        Data(1).Cell_DFoF(frame,cell) = (currentInt-initialInt)/(initialInt);
    end
end
toc;
stopclock = stopclock + toc;
%% Calculate FFT Based on Cell dF/F 
tic;
Data.Cell_FFT = [];

%FFT Parameters
L = length(F_RANGE);
t = (0:L-1)/fs;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = fs/2*linspace(0,1,NFFT/2+1); %x-axis (frequency axis) for plotting


%calculate Delta F of F for Cells
for cell = C_RANGE
    Data.Cell_FFT(:,cell) = fft ( Data.Cell_DFoF(:,cell) , NFFT) / L ;
end
toc;
stopclock = stopclock + toc;
%% Cell Output Range
C_RANGE_OUTPUT = C_RANGE; %[31,101,231,161,171,261]; %(C_RANGE ~= (1) & C_RANGE ~= (1)) .* C_RANGE;

%% Outputs
if 1>0
    %% Output Matrix
    
C_RANGE_OUTPUT(C_RANGE_OUTPUT==0) = []; 

C_RANGE_OUTPUT = unique(C_RANGE_OUTPUT);
%remove empty rows of cells not analyzed 
Data.Cell_DFoF_WithoutEmpt = [];
cnt = 1;
for cell = C_RANGE_OUTPUT
     if sum(Data.Cell_DFoF(:,cell))~=0
         Data.Cell_DFoF_WithoutEmpt(:,cnt) = Data.Cell_DFoF(:,cell);
         cnt = cnt + 1;
     end
end

stopclock_finaltime = duration(seconds(stopclock),'Format','mm:ss');
disp('Total Run Time: '); 
disp(stopclock_finaltime);
    %% previews
    close all;
    fig(1) = figure('Name','Overview of ROI Creation');
subplot(331); imagesc(im); axis 'image'; title("Initial Image");
subplot(332); imagesc(im0); axis 'image'; title("ImAdjusted"); 
subplot(333); imagesc(mask); axis 'image'; title("Separating ROI from BG");
subplot(334); imagesc(im_Filt); axis 'image'; title("Filtered and ImAdjusted"); 
subplot(335); imagesc(im_Skel); axis 'image'; title("Segmentation");
subplot(336); imagesc(im_Inv); axis 'image'; title("Inverted"); 
subplot(337); imagesc(im_Cull); axis 'image'; title("Culled large Areas"); 
subplot(338); imagesc(labl); axis 'image'; title("Cells/Objects Identified: "+num_labl);
subplot(339); imagesc(imdilate(im_Cull, SEdilate) - im_Cull); axis 'image'; title("Dilated Edges"); 

fuseROI = imfuse(imadjust(im),imdilate(im_Cull, SEdilate) - im_Cull);
fuseBG_ROI = imfuse(imadjust(im),Max_BG_ROI);

fig(2) = figure('Name','Overview of ROIs'); 
subplot(121); imagesc(fuseROI); axis 'image'; title("Cell/Obj ROI Preview"); 
subplot(122); imagesc(fuseBG_ROI); axis 'image'; title("Max BG ROI Preview"); 

movegui(fig(1), 'northeast');
movegui(fig(2), 'northwest');

pause(0.00001);
    %% Plot BG, Mean, Indiv, and FFT of dFoF 
movingmean_wndw = 1;  
    
col = jet(length(C_RANGE_OUTPUT));
% set(0,'DefaultAxesColorOrder',col); 
whitebg([0.8 0.8 0.8]);
%set(gca,'NextPlot','replaceall','ColorOrder',copper);

fig(3) = figure('Name','BG,Mean, and All Individual DeltaFoF plots and FFT plots'); movegui(fig(3), 'north')
subplot(221); plot(movmean(Data.BG_DFoF,movingmean_wndw)); title("BG \DeltaF/F"); xlim([F_RANGE(1), max(F_RANGE)]); ylim([-axis_y_lim axis_y_lim]); axis 'square'; xlabel('Frame'); ylabel('\DeltaF/F');
subplot(223); plot(movmean(mean(Data.Cell_DFoF_WithoutEmpt,2),movingmean_wndw)) ; title("Mean \DeltaF/F of "+length(C_RANGE_OUTPUT)+" Selected Cells/Objs"); xlim([F_RANGE(1), max(F_RANGE)]); ylim([-axis_y_lim axis_y_lim]); axis 'square'; xlabel('Frame'); ylabel('\DeltaF/F');

%Indiv Plots on a single
subplot(224);
m_cnt = length(C_RANGE_OUTPUT);
cnt = 1;
while cnt < m_cnt
    hold on
    plot(movmean(Data.Cell_DFoF_WithoutEmpt(:,cnt),movingmean_wndw),line_style,'Color',col(cnt,:)); 
    cnt = cnt + 1;
end
title("Individual \Delta F/F of "+length(C_RANGE_OUTPUT)+" Selected Cells/Objs"); xlim([F_RANGE(1), max(F_RANGE)]); ylim([-axis_y_lim axis_y_lim]); axis 'square'; xlabel('Frame'); ylabel('\DeltaF/F');
hold off
legend(strsplit(num2str(C_RANGE_OUTPUT)));
legend('hide');
datacursormode on;
dcm = datacursormode(fig(3));
set(dcm,'UpdateFcn',@customdatatip);


%FFT
subplot(222);
m_cnt = length(C_RANGE_OUTPUT);
cnt = 1;
while cnt < m_cnt
    hold on
    plot(f,2*abs(Data.Cell_FFT(1:NFFT/2+1,C_RANGE_OUTPUT(cnt))),'Color',col(cnt,:)); 
    cnt = cnt + 1;
end
title("FFT of "+length(C_RANGE_OUTPUT)+" Selected Cells/Objs"); xlabel('Frequency (Hz)'); ylabel('|Y(f)|'); xlim([0.1 fs/2]);
legend(strsplit(num2str(C_RANGE_OUTPUT)));
legend('hide');
hold off
    %% Plot Kymographs
%Kymographs
fig(4) = figure('Name','Kymographs of DeltaFoF'); movegui(fig(4), 'south'); 

%Custom Colormap for the kymograph
hi_col = [0 0 1];               % the positive values
lo_col = [1 0 0];               % the negative values
mid_col = [0.9 0.9 0.9];        % the mid value color

resol = 128;
cmapkymotop = [linspace(mid_col(1),hi_col(1),resol)',linspace(mid_col(2),hi_col(2),resol)',linspace(mid_col(3),hi_col(3),resol)'];
cmapkymobot = [linspace(lo_col(1),mid_col(1),resol)',linspace(lo_col(2),mid_col(2),resol)', linspace(lo_col(3),mid_col(3),resol)'];
cmapKymo = [cmapkymobot;cmapkymotop]; 

resol1 = 120;
s_f = 1;
custom1 = [linspace(0,0,resol1)',linspace(0,1,resol1)',linspace(1,1,resol1)'];
custom2 = [linspace(0.75,0,resol1*s_f)',linspace(0.75,0,resol1*s_f)',linspace(0.75,1,resol1*s_f)'];
custom3 = [linspace(1,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)'];
custom4 = [linspace(1,1,resol1)',linspace(1,0,resol1)',linspace(0,0,resol1)'];
%[custom3;custom2;custom1]
%colormap([custom4;custom3;custom2;custom1]) ; caxis([-1 1]);
cmaprange = [custom4;custom3;custom2;custom1];

%Subplots for the various kymographs to be created 
colormap(cmaprange);
subplot(4,3,[1 3]); imagesc(Data.Cell_DFoF_WithoutEmpt'); axis 'tight'; xlim([F_RANGE(1), max(F_RANGE)]); title("Individual \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel('Cell No.'); yticks(1:length(C_RANGE_OUTPUT)); yticklabels(C_RANGE_OUTPUT); caxis([-1 1]); %colorbar;
subplot(4,3,[4 6]); imagesc(repmat(mean(Data.Cell_DFoF_WithoutEmpt,2)',20,1)); axis 'image'; xlim([F_RANGE(1), max(F_RANGE)]); title("Average \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel("Avg of "+length(C_RANGE_OUTPUT)+" Cells"); set(gca,'TickLength',[0 0],'YTickLabel',[]); caxis([-1 1]); %colorbar; 
subplot(4,3,[7 9]); imagesc(diff(Data.Cell_DFoF_WithoutEmpt)'); axis 'tight'; xlim([F_RANGE(1), max(F_RANGE)]); title("Individual Difference \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel('Cell No.'); yticks(1:length(C_RANGE_OUTPUT)); yticklabels(C_RANGE_OUTPUT); caxis([-1 1]); %colorbar;
subplot(4,3,[10 12]); imagesc(repmat(mean(diff(Data.Cell_DFoF_WithoutEmpt),2)',20,1)); axis 'image'; xlim([F_RANGE(1), max(F_RANGE)]); title("Average Difference \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects");xlabel('Frame'); ylabel("Avg of "+length(C_RANGE_OUTPUT)+" Cells"); set(gca,'TickLength',[0 0],'YTickLabel',[]); caxis([-1 1]); %colorbar; 
axes('Position', [0.085 0.55 0.9 0.4], 'Visible', 'off'); h1 = colorbar; caxis([-1 1]); 
axes('Position', [0.085 0.11 0.9 0.4], 'Visible', 'off'); h2 = colorbar; caxis([-1 1]);
    %% Overview of Cells Selected
    
    %Original Image to overlay with 
    main_image = Storage(F_RANGE(1)).OG_im;
    
    %Create ROI with selected objects/cells
    blnk = zeros(im_height,im_width);
    for cell = C_RANGE_OUTPUT
        add_to = (Storage(F_RANGE(1)).Objs==cell)*cell;
        blnk = blnk + add_to;
    end
    
    %dialte ROI of objects/cells as above
    blnk_di = imdilate(blnk,SEdilate)-blnk;
    
    %overlay dialted ROI and main image
    composite = imfuse(imadjust(main_image), blnk_di>0);
    
    %Output figure
    fig(5) = figure('Name','Overview of Labeled Cell-Objs');
    movegui(fig(5), 'northwest')
    subplot(121); imagesc(composite); axis 'image';
    
    %Number the cells with their proper number
    s = regionprops(logical(Storage(F_RANGE(1)).Objs), 'Centroid');
    hold on
        for k = C_RANGE_OUTPUT
            c = s(k).Centroid;
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'color','cyan');
        end
        hold off;
        
        
%     blnk = zeros(im_height,im_width);
%     for cell = C_RANGE_OUTPUT
%         add_to = Storage(F_RANGE(1)).Objs==cell;
%         blnk = blnk + add_to;
%     end

darkest_val = [0,0,0];
col_f = [darkest_val;col];
       subplot(122); colormap(col_f); imagesc(blnk);axis('image')
    %% Plot Cell Statistics/Metrics
fig(6) = figure('Name','Cell-Obj Statistics'); movegui(fig(6), 'southeast')
subplot(231); hist([stats(C_RANGE_OUTPUT).Area],25); title("Area of "+length(C_RANGE_OUTPUT)+" Cells/Objs"); axis square; xlabel("Area Size(px)");
subplot(232); hist([stats(C_RANGE_OUTPUT).Perimeter],25); title ("Perimeter of "+length(C_RANGE_OUTPUT)+" Cells/Objs"); axis square; xlabel("Perimeter Length(px)");
subplot(233); bar([stats(C_RANGE_OUTPUT).Eccentricity]); title ("Eccentricity of "+length(C_RANGE_OUTPUT)+" Cells/Objs"); axis square; xlabel("Cell Number");
subplot(235); hist( abs(deg2rad([stats(C_RANGE_OUTPUT).Orientation])) ,20 ); title ("Orientation of "+length(C_RANGE_OUTPUT)+" Cells/Objs"); axis square; 
    %% Indiv Cell DFoF Plots 
%Indiv Plots
Figure_on_current_page = -1*indvplot_x*indvplot_y;
        
        page_val = 1;
        
        count_total_graphs = size(Data.Cell_DFoF_WithoutEmpt,2);
        current_graph = 1;
        
        while current_graph < count_total_graphs+1
            if  Figure_on_current_page ~= 0
                fig(6+page_val) = figure('Name',"Individual Cell-Objs on Own Plot - Page "+page_val+" of "+(ceil(count_total_graphs/(indvplot_x*indvplot_y))));
                movegui('southwest')
                while Figure_on_current_page ~= 0
                %plot options and preferences 
                    subplot(indvplot_y,indvplot_x,Figure_on_current_page+(indvplot_x*indvplot_y+1));
                    plot(movmean(Data(1).Cell_DFoF(:,C_RANGE_OUTPUT(current_graph)),movingmean_wndw),'color',col(current_graph,:));
                
                    
                    title ("Cell No. "+C_RANGE_OUTPUT(current_graph));
                    grid on; %grid minor;
                    axis 'square';
                    %xticks(1:floor(max_no_frame/total_video_time_min):max_no_frame);
                    %mean_plot = mean(Vals_Storage(1).Cell_DeltaFoF(:,current_graph));
                    ylim([-axis_y_lim , axis_y_lim]); %y-axis limits
                    xlim([F_RANGE(:,1),max(F_RANGE)]); %x-axis limits
                    xlabel('Frame'); % x-axis label
                    ylabel('\DeltaF/F'); % y-axis label
                %Values to change for Page process
                    Figure_on_current_page = Figure_on_current_page +1;
                    current_graph = current_graph + 1 ;
                    
                    if current_graph > count_total_graphs %stops count from going over
                        break
                    end
                end
                Figure_on_current_page = -1*indvplot_x*indvplot_y;
                page_val = page_val +1;
            end
        end 
end
%%

return


%% Save Figures

%create Folder names depending on cell output
if C_Analysis_Percent == 100 %&& C_RANGE_OUTPUT == C_RANGE
    cells_name = char("All_Cells_from+" + fileX);
elseif C_Analysis_Percent < 100 %&& C_RANGE_OUTPUT == C_RANGE
    cells_name = char(C_Analysis_Percent + "%_of_All_Cells_from_"+fileX);
end
destFolder = fullfile('Matlab_Analysis/' , cells_name);

%make these folders
mkdir(pathX, destFolder);

FolderName = fullfile(pathX,destFolder);   % Your destination folder

%save raw data into folder
    Analysis_Paramenters.C_Analysis_Percent = C_Analysis_Percent;
    Analysis_Paramenters.dilate_size = dilate_size;
    Analysis_Paramenters.Thresh_Mod = Thresh_Mod;
    Analysis_Paramenters.Clean_Mod = Clean_Mod; 
    Analysis_Paramenters.Purne_Extensions_Mod = Purne_Extensions_Mod;
    Analysis_Paramenters.Min_Obj_Size = Min_Obj_Size;
    Analysis_Paramenters.Max_Obj_Size = Max_Obj_Size;
    Analysis_Paramenters.F_RANGE = F_RANGE;
    Analysis_Paramenters.C_RANGE = C_RANGE;
    Analysis_Paramenters.C_RANGE_OUTPUT = C_RANGE_OUTPUT;
    Analysis_Paramenters.C_RANGE = C_RANGE;
    Analysis_Paramenters.C_Analysis_Percent = C_Analysis_Percent;
    Analysis_Paramenters.NFFT = NFFT;
    Analysis_Paramenters.fs = fs; 
    Analysis_Paramenters.axis_y_lim = axis_y_lim;
    Analysis_Paramenters.indvplot_x = indvplot_x;
    Analysis_Paramenters.indvplot_y = indvplot_y;
    
    save( fullfile(FolderName, char("Data From-"+fileX+".mat")) , 'Data') ;
    save( fullfile(FolderName, char("Analysis Parameters From-"+fileX+".mat")), 'Analysis_Paramenters');

%save figures into folders
for iFig = 1:length(fig)
  FigName = fig(iFig).Name;
  savefig(fig(iFig), fullfile(FolderName, FigName));
  disp("Completed Saving Figure " +iFig); 
end

%%
Storage(1).OG_im = imread('/Volumes/KOS-Extra-1/Losert Lab/Microscope Raw Data/Multiscale /2018_08_22/HEK High Speed 2018_08_22/1x1FullField_70.5fps_3/1x1FullField_70.5fps_3_MMStack_Pos0.ome.tif',1);

%% CLEAR and CLOSE ALL
clear; close all;
%% play movie

for f = F_RANGE
    figure(36538);colormap(copper); imagesc(imgaussfilt((Storage(f).OG_im),1)); title("Frame:"+f);
    pause(0.2)
end
%% play movie
for f = F_RANGE
    figure(36538);colormap(copper); imagesc(imgaussfilt((imread(im_path,f)),25)); colorbar 
    pause(0.1)
end


%% create diff image to see waves
for frame = F_RANGE(2):(length(F_RANGE)-1)
    diffim = imsubtract( Storage(frame).OG_im , Storage(frame-1).OG_im);
    %Storage(frame).Diff = diffim;
    Storage(frame).Final = imfuse(Storage(frame).OG_im, diffim);
end

%%
for f = F_RANGE(2):(length(F_RANGE)-1)
minsf(f) = min(Storage(f).Final(:));
maxsf(f) = max(Storage(f).Final(:));
end;

min_valf = min(minsf);
max_valf = max(maxsf);

%% create dF/F image
for frame = F_RANGE(2):(length(F_RANGE)-1)
    dFoFim = ((Storage(frame).OG_im) - Storage(F_RANGE(1)).OG_im) ;%./ (Storage(F_RANGE(1)).OG_im);
    Storage(frame).ImdFoF = dFoFim; %imfuse(Storage(frame).OG_im, diffim);
    
end
%% play diff image movie
for f = F_RANGE(2):(length(F_RANGE)-1)
    figure(36538); imagesc(Storage(f).Final); title("Frame "+f); %caxis([min_valf max_valf]);
    pause(0.1)
end

%% save video
frm_rt = 10;
v = VideoWriter(char("video_from_"+fileX+"_at_"+frm_rt+"_fps.avi"));
v.FrameRate = frm_rt;
open(v);
for frame = F_RANGE(2):(length(F_RANGE)-1)
   frame_v = Storage(frame).Final;
   writeVideo(v,frame_v);
   
end

close(v);

%You must call open(v) before calling writeVideo.

%% 
figure;
for frame = F_RANGE(2):(length(F_RANGE)-1)
    imagesc(Storage(frame).Diff); axis 'image'; title("frame"+frame); caxis([min_valf max_valf]);
    pause(0.1);
end
%%

imagesc(Storage(24).ImDiff>250);
%%
spy(sparsce(Storage(24).ImDiff>250));
%%

for f = F_RANGE(2):(length(F_RANGE)-1)

    mm = sparse(Storage(f).ImDiff>220); % save the binarized image as sparse matrix with coordinates for ever 1
    [ xcoor,ycoor,valu ] = find(mm); %store x and y coordinates and values
    p = polyfit(xcoor,ycoor,3);
    %x1 = linspace(0,5);
    %y1 = linspace(0,nImage);
    p1 = polyval(p,1:max(ycoor));
    spy(sparse(mm)); %lsline;
   
    hold on; plot(p1);
    pause(0.2); %hold off;
end 

%%

 %% Plot Kymographs
%Kymographs
fig(4) = figure('Name','Kymographs of DeltaFoF'); movegui(fig(4), 'south'); 

%Custom Colormap for the kymograph
hi_col = [0 0 1];               % the positive values
lo_col = [1 0 0];               % the negative values
mid_col = [0.9 0.9 0.9];        % the mid value color

resol = 128;
cmapkymotop = [linspace(mid_col(1),hi_col(1),resol)',linspace(mid_col(2),hi_col(2),resol)',linspace(mid_col(3),hi_col(3),resol)'];
cmapkymobot = [linspace(lo_col(1),mid_col(1),resol)',linspace(lo_col(2),mid_col(2),resol)', linspace(lo_col(3),mid_col(3),resol)'];
cmapKymo = [cmapkymobot;cmapkymotop]; 

resol1 = 120;
s_f = 1;
custom1 = [linspace(0,0,resol1)',linspace(0,1,resol1)',linspace(1,1,resol1)'];
custom2 = [linspace(0.75,0,resol1*s_f)',linspace(0.75,0,resol1*s_f)',linspace(0.75,1,resol1*s_f)'];
custom3 = [linspace(1,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)',linspace(0,0.75,resol1*s_f)'];
custom4 = [linspace(1,1,resol1)',linspace(1,0,resol1)',linspace(0,0,resol1)'];
%[custom3;custom2;custom1]
%colormap([custom4;custom3;custom2;custom1]) ; caxis([-1 1]);
cmaprange = [custom4;custom3;custom2;custom1];

%Subplots for the various kymographs to be created 
colormap(cmaprange);
subplot(4,3,[1 3]); imagesc(Data.Cell_DFoF_WithoutEmpt'); axis 'tight'; xlim([F_RANGE(1), max(F_RANGE)]); title("Individual \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel('Cell No.'); yticks(1:length(C_RANGE_OUTPUT)); yticklabels(C_RANGE_OUTPUT); caxis([-1 1]); %colorbar;
subplot(4,3,[4 6]); imagesc(repmat(mean(Data.Cell_DFoF_WithoutEmpt,2)',20,1)); axis 'image'; xlim([F_RANGE(1), max(F_RANGE)]); title("Average \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel("Avg of "+length(C_RANGE_OUTPUT)+" Cells"); set(gca,'TickLength',[0 0],'YTickLabel',[]); caxis([-1 1]); %colorbar; 
subplot(4,3,[7 9]); imagesc(diff(Data.Cell_DFoF_WithoutEmpt)'); axis 'tight'; xlim([F_RANGE(1), max(F_RANGE)]); title("Individual Difference \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects"); xlabel('Frame'); ylabel('Cell No.'); yticks(1:length(C_RANGE_OUTPUT)); yticklabels(C_RANGE_OUTPUT); caxis([-1 1]); %colorbar;
subplot(4,3,[10 12]); imagesc(repmat(mean(diff(Data.Cell_DFoF_WithoutEmpt),2)',20,1)); axis 'image'; xlim([F_RANGE(1), max(F_RANGE)]); title("Average Difference \DeltaF/F of "+size(C_RANGE,2)+" Cells/Objects");xlabel('Frame'); ylabel("Avg of "+length(C_RANGE_OUTPUT)+" Cells"); set(gca,'TickLength',[0 0],'YTickLabel',[]); caxis([-1 1]); %colorbar; 
axes('Position', [0.085 0.55 0.9 0.4], 'Visible', 'off'); h1 = colorbar; caxis([-1 1]); 
axes('Position', [0.085 0.11 0.9 0.4], 'Visible', 'off'); h2 = colorbar; caxis([-1 1]);
%% CUSTOM FUNCTION - ALLOWS FOR Data Cursor update
function output_txt = customdatatip(obj,event_obj,str)
pos = get(event_obj, 'Position');
output_txt = {...
    ['X: ', num2str(pos(1),4)]...
    ['Y: ', num2str(pos(2),4)] ...
    ['Cell/Obj No: ', event_obj.Target.DisplayName]...
};
% from https://stackoverflow.com/questions/31244981/show-matlab-legend-text-in-datatip?rq=1
end
