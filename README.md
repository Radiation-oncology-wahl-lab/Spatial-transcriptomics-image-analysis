# Spatial-transcriptomics-image-analysis
Spatial transcriptomics Immunofluorescence(IF) image analysis and annotation (compatible with 10x Genomics visium and Cytassist).
This code was developed on Matlab and uses the Image Processing Toolbox functions.
This code is meant for use with conventional 10x Genomics visium and Cytassist slides after immunofurescent staining. However, it may be user defined and can be expanded to include analysis of IF images paired with Visium and Cytassist HD. The analysis output can be imported to Loupe Browser or be used with other bioinformatic analysis tools and programs.

The code allows the user to objectively annotate the barcodes (spots in case of convetional visium and cytassist sldes) based on fluorescence intensity and percentage positivity of the chosen biomarker.

Create the following 5 scripts in matlab:

1) Cell_Sweep.m:
function [Cell_table, non_binary_img, row, Cell_pos, Numb_total] = Cell_Sweep(non_binary_img, red_channel, blue_channel, green_channel, Perimeter_cell, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup)
      
   Size_non_binary = size(non_binary_img); 
   
   %Making new image gray, binary, getting rid of noise, etc. 
    Green_gray_cell_img = rgb2gray(non_binary_img);
    Green_binary_cells_img = imbinarize(Green_gray_cell_img);
    Pixel_Removal_Less_Twenty = bwareaopen(Green_binary_cells_img, Pixel_cleanup); 
    
    %Getting cell data 
    Cell_regions = regionprops(Pixel_Removal_Less_Twenty, 'centroid', 'Perimeter', 'PixelList', 'MajorAxisLength');
    
    Size_regions = size(Cell_regions);
    
    if Size_regions(1) == 0 
        Cell_table = 0; 
        row = 0;
        Cell_pos = 0; 
        Numb_total = 0; 
        return; 
    end 
    
    Perimeter = cat(2, Cell_regions.Perimeter); 
    Number_tot_cells = size(Perimeter); 
    Numb_total = Number_tot_cells(2); 
    
    row = 0;
    Cell_pos = 0; 
    
    Cell_table(1000) = struct('centroids', [], 'positivity', []);
    
    %Getting info of cells that meet cell requirements
    for indx_solid = 1: Number_tot_cells(2)
        if Perimeter(indx_solid) <= Perimeter_cell  
            row = row + 1;  
            
            Ex_Pixels = cat(3, Cell_regions(indx_solid).PixelList);  
            
            Size_pixels = size(Ex_Pixels);
            
            alpha = alphaShape(Ex_Pixels(:,1), Ex_Pixels(:,2));
            [~ , boundary_coord] = boundaryFacets(alpha); 
             
            binary_poly = roipoly([1 Size_non_binary(1)], [1 Size_non_binary(2)], non_binary_img, boundary_coord(:,1), boundary_coord(:,2)); 
            
            boundary_mask = boundarymask(binary_poly, 4); 
            cell_mask = imfill(boundary_mask, 'holes'); 
            
            Mask_matrix = cat(3, cell_mask, cell_mask, cell_mask); 
            non_binary_img(Mask_matrix) = 0;
            
            Sum_positive = sum(sum(red_channel(binary_poly) >= Red_cutoff & blue_channel(binary_poly) >= Blue_cutoff & green_channel(binary_poly) >= Green_cutoff)); 
            Percent_pos = (Sum_positive/(Size_pixels(1))) * 100; 
            
                if Percent_pos >= Percent_cutoff  
                    Cell_table(row).positivity = 'pos';
                    Cell_pos = Cell_pos + 1; 
                else 
                    Cell_table(row).positivity = 'neg'; 
                end 

            Cell_table(row).centroids = Cell_regions(indx_solid).Centroid;

        end 
    end
end

2) Cluster_array.m:
function [cluster_images, Circle_matrix, Tissue_Img, Final_Cluster_table, Zero_table, Cell_table] = Cluster_array(x, y, Tissue_Img, Circle_matrix, radius, range, Perimeter_cell, Clusters_orig, Program_call, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup) 

tic;

Circle_matrix(:,3) = Circle_matrix(:,3) + 11; 

%Creating clusters with given radius
Predata = linspace(0, 2*pi, 10000); 
    for indx = 1:size(Circle_matrix,1)
        disp("Cluster" + indx); 
        x_circle_eqn = cos(Predata) * radius + x(indx);
        y_circle_eqn = sin(Predata) * radius + y(indx);
        Clusters = poly2mask(x_circle_eqn, y_circle_eqn, size(Tissue_Img,1), size(Tissue_Img,2));
        mask = bsxfun(@times, Tissue_Img, cast(Clusters, class(Tissue_Img)));
        props = regionprops(Clusters, 'BoundingBox');
        mask = imcrop(mask, props.BoundingBox);
        cluster_images{indx} = mask;
    end
    
    Cluster_headings = {'x coordinate', 'y coordinate', 'Number Cells', 'Number Positive', 'Intensity', 'Percent Positive', 'Program Call'};
    Final_Cluster_table = zeros(size(Circle_matrix,1) + 1, 7);
    
    for indx_cluster = 1:size(Circle_matrix,1)
        disp("Cluster" + indx_cluster); 
        [Total_pos, Total_cells, Red_intensity, Cell_table] = Cluster_Pos(cluster_images{indx_cluster}, size(Circle_matrix, 1), range, radius, Perimeter_cell, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup); 
        Final_Cluster_table(indx_cluster, 1) = x(indx_cluster); 
        Final_Cluster_table(indx_cluster, 2) = y(indx_cluster); 
        Final_Cluster_table(indx_cluster, 3) = Total_cells; 
        Final_Cluster_table(indx_cluster, 4) = Total_pos; 
        Final_Cluster_table(indx_cluster, 5) = Red_intensity; 
        Final_Cluster_table(indx_cluster, 6) = ((Total_pos)/(Total_cells)) * 100; 
        if Program_call == 'a'
            if Red_intensity >= Intensity_cutoff 
                Final_Cluster_table(indx_cluster, 7) = 1;
            else 
                Final_Cluster_table(indx_cluster, 7) = 0; 
            end 
        end 
        if Program_call == 'b'
            if (((Total_pos)/(Total_cells)) * 100) >= Percent_cutoff
                Final_Cluster_table(indx_cluster, 7) = 1;
            else
                Final_Cluster_table(indx_cluster, 7) = 0; 
            end 
        end 
        if Program_call == 'c'
            if Red_intensity >= Intensity_cutoff && (((Total_pos)/(Total_cells)) * 100) >= Percent_cutoff
                Final_Cluster_table(indx_cluster, 7) = 1;
            else
                Final_Cluster_table(indx_cluster, 7) = 0;
            end 
        end 
        
    end
    
    Zero_table = zeros(size(Circle_matrix,1), 7);
    One_table = zeros(size(Circle_matrix,1), 7); 
    Two_table = zeros(size(Circle_matrix,1), 7); 
    Three_table = zeros(size(Circle_matrix,1), 7); 
    Four_table = zeros(size(Circle_matrix,1), 7); 
    Error_table = zeros(size(Circle_matrix,1), 7); 
    
    
    zero = 0;
    one = 0;
    two = 0;
    three = 0; 
    four = 0; 
    error = 0; 
    for indx_cluster = 1:size(Circle_matrix,1)   
        if Final_Cluster_table(indx_cluster, 7) == 0
            zero = zero + 1; 
            Zero_table(zero, :) = Final_Cluster_table(indx_cluster, :); 
        elseif Final_Cluster_table(indx_cluster, 7) == 1
            one = one + 1; 
            One_table(one, :) = Final_Cluster_table(indx_cluster, :); 
        elseif Final_Cluster_table(indx_cluster, 7) == 2
            two = two + 1; 
            Two_table(two, :) = Final_Cluster_table(indx_cluster, :);
        elseif Final_Cluster_table(indx_cluster, 7) == 3
            three = three + 1; 
            Three_table(three, :) = Final_Cluster_table(indx_cluster, :);
        elseif Final_Cluster_table(indx_cluster, 7) == 4
            four = four + 1; 
            Four_table(four, :) = Final_Cluster_table(indx_cluster, :);
        else
            error = error + 1; 
            Error_table(error,:) = Final_Cluster_table(indx_cluster, :); 
        
        end
    end 
    
    Final_Cluster_table = [Cluster_headings; num2cell(Final_Cluster_table)]; 
    
    imshow(Clusters_orig);
    hold on; 
    plot(One_table(:,1), One_table(:,2), 'm*');
    hold off; 

    toc; 
end

3) Cluster_pos.m:
function [Total_pos, Total_cells, Red_intensity, Cell_table] = Cluster_Pos(crop, Cluster_Num, range, radius, Perimeter_cell, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup)

%Getting Data for Crop
Int_Vals = crop; 

%Color Channels
Red = Int_Vals(:,:,1);
Blue = Int_Vals(:,:,3);
Green = Int_Vals(:,:,2); 

if Intensity_color == "Red" 
    Intensity_color = Red; 
elseif Intensity_color == "Blue"
    Intensity_color = Blue;
elseif Intensity_color == "Green"
    Intensity_color = Green;
end

%Preallocating Data
size_struct = 1500; 
All_Cell_table(size_struct) = struct('centroids', []);
adjust = 0;  
Cluster_matrix = zeros(Cluster_Num, 3); 

crop_size = size(Int_Vals);
x_center = round(crop_size(1)/2);
y_center = round(crop_size(2)/2);

Green_cells_img = Int_Vals;  
Total_pos = 0; 
Total_cells = 0; 
Red_intensity = 0; 
    
    for int_pixel_thresh = range

        Blue_thresh = Blue >= int_pixel_thresh; 
        Blue_mask = cat(3, Blue_thresh, Blue_thresh, Blue_thresh); 
        Green_cells_img(~Blue_mask) = 0; 
        
        if int_pixel_thresh == range(1) 
                Red_cells_img = Green_cells_img; 
                Red_thresh = Intensity_color >= Intensity_cutoff; 
                Red_mask = cat(3, Red_thresh, Red_thresh, Red_thresh); 
                Red_cells_img(~Red_mask) = 0; 
                Red_final = Red_cells_img(:,:,1); 
                Sum_cluster_red_intensity = sum(sum(Red_final)); 
                Total_cluster_pixels = sum(sum(Red_final >= Intensity_cutoff)); 
                Red_intensity = Sum_cluster_red_intensity/Total_cluster_pixels; 
                if Total_cluster_pixels == 0
                    Red_intensity = 0;
                end
        end
    
        [Cell_table, img, row, Numb_pos, Cell_numb] = Cell_Sweep(Green_cells_img, Red, Blue, Green, Perimeter_cell, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup); 
    
        Green_cells_img = img; 
    
        ID_Cell = Cell_table;
        
        Total_pos = Total_pos + Numb_pos; 
        Total_cells = Total_cells + Cell_numb; 
    
        for indx_row = 1:row 
            All_Cell_table(indx_row + adjust).centroids = ID_Cell(indx_row).centroids; 
            All_Cell_table(indx_row + adjust).positivity = ID_Cell(indx_row).positivity;
        end 
       
        adjust = adjust + row; 
    end 
            
        
 
end 
    
4) Cluster_Radius.m:
function [Clusters_orig, Circle_matrix, range, radius, Perimeter_cell, Program_call, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup] = Cluster_Radius(x, y, Tissue_Img)

%User inputs 
radius = 99;
range = []; 
Cluster_number = 0; 
Perimeter_cell = []; 
Program_call = [];
Intensity_color = []; 
Intensity_cutoff = [];
Percent_cutoff = [];
Red_cutoff = [];
Blue_cutoff = [];
Green_cutoff = []; 
Pixel_cleanup = [];
dlgtitle_2 = 'Input for Cluster Positivity Program';
prompt_2 = {'Please enter the range for cell sweeping:', 'Please enter the number of clusters:', 'Average Cell Perimeter:', 'Please Enter Program Call Choice (Intensity (a), Percent (b), or Both (c)):', 'Intensity Color (Red, Green, or Blue)(Leave as Default if (b) Chosen):', 'Intensity Lower Threshold (Leave as Default if (b) Chosen):', 'Percent Positive Threshold (Leave as Default if (a) Chosen):', 'Red Lower Pixel Value Threshold (Leave as Default if (a) Chosen):', 'Blue Lower Pixel Value Threshold (Leave as Default if (a) Chosen):', 'Green Lower Pixel Value Threshold (Leave as Default if (a) Chosen):', 'Number of Pixels to not be Considered a Cell:'};
dialog_size_2 = [1 75];
default_input_2 = {'50:25:250', '0', '80', 'a', 'Red', '55', '5', '55', '50', '0', '20'};

    while Cluster_number == 0
        input_info = inputdlg(prompt_2, dlgtitle_2, dialog_size_2, default_input_2);
        range = str2num(input_info{1});  
        Cluster_number = str2double(input_info{2});
        Perimeter_cell = str2double(input_info{3});
        Program_call = input_info{4};
        Intensity_color = input_info{5};
        Intensity_cutoff = str2num(input_info{6});
        Percent_cutoff = str2num(input_info{7});
        Red_cutoff = str2num(input_info{8});
        Blue_cutoff = str2num(input_info{9});
        Green_cutoff = str2num(input_info{10});
        Pixel_cleanup = str2num(input_info{11});

        if (isempty(input_info))
            uiwait(msgbox('All values must be entered to run the program'));
            continue;
        end
    end
    
    Radius_matrix = zeros(size(x,1), 1); 
    Radius_matrix(:,1) = radius; 
    Circle_matrix = [x, y, Radius_matrix];

    %Creating clusters with given radius
    Clusters_orig = insertShape(Tissue_Img, "FilledCircle", Circle_matrix, 'Color', [255 255 255], 'Opacity', 0.2); 


end


5) Cluster_Positivity_Run.m:
%Cluster Positivity 
 %Based on specifications by user, program places cell clusters on tif
 %image and finds/analyzes cells to determine cell positivity for
 %Gamma-H2AX

Radius_check = []; 
Cluster_check = 0; 
Whole_Image_Name = [];
Excel_File_Name = [];
Img_Name = '___________.tif';
Excel_Name = '___________.csv';

dlgtitle = 'Input for File Names';
prompt = {'Please enter in the image filename:','Please enter Excel file name containing cluster coordinates:'};
dialog_size = [1 75];
default_input = {Img_Name, Excel_Name};

input_info = inputdlg(prompt, dlgtitle, dialog_size, default_input);

while isempty(input_info) 
    if (isempty(input_info))
            uiwait(msgbox('All values must be entered to run the program'));
            continue;
    end
end 


Whole_Image_Name = input_info{1}; 
Excel_File_Name = input_info{2}; 

%Reading in image and Excel file containing cluster coordinates
Tissue_Img = imread(Whole_Image_Name); 
coordinates = csvread(Excel_File_Name, 1, 1);

x = coordinates(:,1); 
y = coordinates(:,2); 

[Clusters_orig, Circle_matrix, range, radius, Perimeter_cell,Program_call, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup] = Cluster_Radius(x, y, Tissue_Img);
    


[cluster_images, Circle_matrix, Tissue_Img, Final_Cluster_table, Zero_table] = Cluster_array(x, y, Tissue_Img, Circle_matrix, radius, range, Perimeter_cell, Clusters_orig, Program_call, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup);

Cluster_array:
function [cluster_images, Circle_matrix, Tissue_Img, Final_Cluster_table, Zero_table, Cell_table] = Cluster_array(x, y, Tissue_Img, Circle_matrix, radius, range, Perimeter_cell, Clusters_orig, Program_call, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup) 

tic;

Circle_matrix(:,3) = Circle_matrix(:,3) + 11; 

%Creating clusters with given radius
Predata = linspace(0, 2*pi, 10000); 
    for indx = 1:size(Circle_matrix,1)
        disp("Cluster" + indx); 
        x_circle_eqn = cos(Predata) * radius + x(indx);
        y_circle_eqn = sin(Predata) * radius + y(indx);
        Clusters = poly2mask(x_circle_eqn, y_circle_eqn, size(Tissue_Img,1), size(Tissue_Img,2));
        mask = bsxfun(@times, Tissue_Img, cast(Clusters, class(Tissue_Img)));
        props = regionprops(Clusters, 'BoundingBox');
        mask = imcrop(mask, props.BoundingBox);
        cluster_images{indx} = mask;
    end
    
    Cluster_headings = {'x coordinate', 'y coordinate', 'Number Cells', 'Number Positive', 'Intensity', 'Percent Positive', 'Program Call'};
    Final_Cluster_table = zeros(size(Circle_matrix,1) + 1, 7);
    
    for indx_cluster = 1:size(Circle_matrix,1)
        disp("Cluster" + indx_cluster); 
        [Total_pos, Total_cells, Red_intensity, Cell_table] = Cluster_Pos(cluster_images{indx_cluster}, size(Circle_matrix, 1), range, radius, Perimeter_cell, Intensity_color, Intensity_cutoff, Percent_cutoff, Red_cutoff, Blue_cutoff, Green_cutoff, Pixel_cleanup); 
        Final_Cluster_table(indx_cluster, 1) = x(indx_cluster); 
        Final_Cluster_table(indx_cluster, 2) = y(indx_cluster); 
        Final_Cluster_table(indx_cluster, 3) = Total_cells; 
        Final_Cluster_table(indx_cluster, 4) = Total_pos; 
        Final_Cluster_table(indx_cluster, 5) = Red_intensity; 
        Final_Cluster_table(indx_cluster, 6) = ((Total_pos)/(Total_cells)) * 100; 
        if Program_call == 'a'
            if Red_intensity >= Intensity_cutoff 
                Final_Cluster_table(indx_cluster, 7) = 1;
            else 
                Final_Cluster_table(indx_cluster, 7) = 0; 
            end 
        end 
        if Program_call == 'b'
            if (((Total_pos)/(Total_cells)) * 100) >= Percent_cutoff
                Final_Cluster_table(indx_cluster, 7) = 1;
            else
                Final_Cluster_table(indx_cluster, 7) = 0; 
            end 
        end 
        if Program_call == 'c'
            if Red_intensity >= Intensity_cutoff && (((Total_pos)/(Total_cells)) * 100) >= Percent_cutoff
                Final_Cluster_table(indx_cluster, 7) = 1;
            else
                Final_Cluster_table(indx_cluster, 7) = 0;
            end 
        end 
        
    end
    
    Zero_table = zeros(size(Circle_matrix,1), 7);
    One_table = zeros(size(Circle_matrix,1), 7); 
    Two_table = zeros(size(Circle_matrix,1), 7); 
    Three_table = zeros(size(Circle_matrix,1), 7); 
    Four_table = zeros(size(Circle_matrix,1), 7); 
    Error_table = zeros(size(Circle_matrix,1), 7); 
    
    
    zero = 0;
    one = 0;
    two = 0;
    three = 0; 
    four = 0; 
    error = 0; 
    for indx_cluster = 1:size(Circle_matrix,1)   
        if Final_Cluster_table(indx_cluster, 7) == 0
            zero = zero + 1; 
            Zero_table(zero, :) = Final_Cluster_table(indx_cluster, :); 
        elseif Final_Cluster_table(indx_cluster, 7) == 1
            one = one + 1; 
            One_table(one, :) = Final_Cluster_table(indx_cluster, :); 
        elseif Final_Cluster_table(indx_cluster, 7) == 2
            two = two + 1; 
            Two_table(two, :) = Final_Cluster_table(indx_cluster, :);
        elseif Final_Cluster_table(indx_cluster, 7) == 3
            three = three + 1; 
            Three_table(three, :) = Final_Cluster_table(indx_cluster, :);
        elseif Final_Cluster_table(indx_cluster, 7) == 4
            four = four + 1; 
            Four_table(four, :) = Final_Cluster_table(indx_cluster, :);
        else
            error = error + 1; 
            Error_table(error,:) = Final_Cluster_table(indx_cluster, :); 
        
        end
    end 
    
    Final_Cluster_table = [Cluster_headings; num2cell(Final_Cluster_table)]; 
    
    imshow(Clusters_orig);
    hold on; 
    plot(One_table(:,1), One_table(:,2), 'm*');
    hold off; 

    toc; 
end

Create a folder with the above scripts along with the images to be processed and their corresponding .csv files containing X,Y coordinates of the visium/Cytassist spots or barcodes. The coordinates file can be obtained from Loupe browser by clicking "Export Projection" listed under "spatial" aspect of Space Ranger analysis.

Set this folder as the working directory in matlab.

Open Cluster_Positivity_Run.m script, chose "Editor" in the tool bar and click "Run".

The user will be prompted to enter the image and the coordinate file names. Be sure to enter the exact name of the image and coordinates file.

Once the files are entered, the users are given option to set the following determinants for the analysis.

Setting	Default Values
(a) Range for Cell Sweeping	50:25:250
  This detects the borders of the cells. The default setting is set at 50: 25: 250 based on prior images analyzed. 
  These values pertain to the DAPI intensity in the image. "50" is the lower limit of DAPI intensity to detect the boundary of a nucleus. "25" is a step increment in DAPI intensity for determining the boundary of the nucleus. This value determines how many iterations are carried out to determine the cell boundaries. "250" is the upper limit of DAPI intenisty. Changing these values will alter the number of cells detected in the spot. 

(b) Number of clusters
  This is the same as the number of barcodes generated for spatial transcriptomic analysis.
  
(c) Average Cell Perimeter	80
  This value was set based on the images analyzed while generating the code. The size of the cells depends on the tissue being studied. This is the number of pixels in a single nucleus. Any DAPI value above the threshold intesity detected outside the set perimeter will be considered as an overlapping nucleus.
  
(d) Intensity Lower Threshold	55
  This pertains to the threshold for detecting presence of the biomarker of interest.

(e) Percent Positive Threshold	5
  This pertains to the minimum percentage of cells within the cluster positive for the biomarker to call the cluster as positive.
  
(f) Red Pixel Value (Lower Threshold) 55
  This pertains to the threshold for detecting presence of the biomarker of interest. If the biomarker shows red fluorescence, this value will be the same as (d). If the marker of interest is yellow, the red value will be lower and the green value will be higher.

(g) Blue Pixel Value (Lower Threshold)	50
  This pertains to the threshold for detecting presence of the biomarker of interest. Since a nuclear marker was used to create this code, underlying blue (DAPI) is present in the image and is hence, set at this value.  
  
(h) Green Pixel Value (Lower Threshold) 0
  This pertains to the threshold for detecting presence of the biomarker of interest. If the biomarker shows red fluorescence, this value will be zero. If the marker of interest is yellow, the red value will be lower and the green value will be higher.
  
(i) Pixel Clean-up	20
  This is the minimum number of pixels to define a nuclues. 

The output is called "Final Cluster table" in the Matlab workspace. This table shows the coordinates of the barcode, total cells in the spot, intensity of the biomarker in the spot, percentage of cells positive for the biomarker in the spot and the final program call deciding whether the spot is positive or negative. This information can be converted to a csv file containing the barcode and annotation as positive or negative. This file can be imported to loupe browser and differential gene expression can be studied. The annotation file may also be used for advanced bioinformatic analysis of the RNA sequencing data in R and/or Python. 

If this code needs to be applied to visium HD, alter the script in Cluster_Radius.m by changing the value of "radius=" to match the pixel size of the RNAseq platform.



