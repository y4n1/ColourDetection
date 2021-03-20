 function result = colourMatrix(filename)    
    rows = 4;
    cols = 4;
    img = imread(filename);  % reading image file
    figure('Name', filename); % display image file
    subplot(2,2,2);
    imshow(img);
    title('Original image');
    img_gray = rgb2gray(img);  % Transform image into BW format
    sz_exp = [480 480 3]; % expected warping size
    [box_size, movingPoint] = box(img);  % Preliminary step
    
    %define moving point & fixed point
    fixedPoint = [26.5,26.5; 445.5,26.5; 445.5,445.5; 26.5,445.5];
    % Swap the order of movingPoint to fit the fixedPoint sequence 
    mov21 = movingPoint(2,1); mov22 = movingPoint(2,2);
    movingPoint(2,1) = movingPoint(3,1); movingPoint(2,2) = movingPoint(3,2);
    movingPoint(3,1) = movingPoint(4,1); movingPoint(3,2) = movingPoint(4,2);
    movingPoint(4,1) = mov21; movingPoint(4,2) = mov22;
    
    % Fit the transformation
    tform = fitgeotrans(movingPoint,fixedPoint,'projective');
    % transform the image with reference 2-D image to world  coordinates
    tform.T
    corrected = imwarp(img,tform,'OutputView',imref2d(sz_exp))
    % --- cropping the corrected image ---
    % define cropping area
    binaryImage_cor = bwareafilt(im2bw(corrected,0.7), 1);
    binaryImage_cor = medfilt2(binaryImage_cor, [1 1]);
    %get all other than black objects
    labeledImage_cor = bwlabel(~binaryImage_cor,8);
    measurements_cor = regionprops(labeledImage_cor, 'BoundingBox', 'MajorAxisLength')
    % Get the biggest size of non black object
    n_max = 1;
    n_start = measurements_cor(1).MajorAxisLength;
    for q = 2:length(measurements_cor)
        if n_start < length(img_gray)
            if (n_start < measurements_cor(q).MajorAxisLength)
                n_max = q;
                n_start = measurements_cor(q).MajorAxisLength;
            end
        else
            n_start = measurements_cor(q).MajorAxisLength;
            n_max = q;
        end    
    end
    % crop the target area
    area_target = measurements_cor(n_max).BoundingBox;
    crop = imcrop(corrected,area_target);
    
    % Assign crop image as process image
    img_process = crop;

    % --- Execute the main function ---
    result = main(img_process, rows,cols);
    
    % Display image
    subplot(2,2,1);
    imshow(img_process);
    title('Colour Target');
    subplot(2,2,3);
    %C = mat2cell(result, [1 1 1 1]);
    uitable('Data',result, 'Position',[5 5 270 200]);
    
    % --- Main function ---
    function matrix_result = main(image, rows, cols)
        count = rows * cols;
        img_gray_target = rgb2gray(image);
        [box_size_target, movingPoint_target] = box(image);
        seg_width = box_size_target(3)/cols;
        seg_height = box_size_target(4)/rows;
        x_start = box_size_target(1);
        hsvImage = rgb2hsv(double(image));
        h = hsvImage(:,:,1);
        s = hsvImage(:,:,2);
        v = hsvImage(:,:,3);
        hsvImage1 = round(hsvImage(:,:,1)*360); % multiply by 360 to get the degree
    
        % subtract red colour 
        red = medfilt2(imsubtract(image(:,:,1), img_gray_target), [2,2]);
        %red = ((hsvImage1<=30)|(hsvImage1>330))& 0.17
        red1 = im2bw(red,0.2); 
        %red1 = medfilt2(red1, [2,2]);
        binaryImage_R = bwareafilt(red1, count);
        labeledImage_R = bwlabel(binaryImage_R);
        % get information of the box [left, top, width, height]
        measurements_R = regionprops(labeledImage_R, 'BoundingBox');
        Pos_R = position(measurements_R, 'r', seg_width, box_size_target(1));
      
        % subtract green colour using HSV
        green1 = ((hsvImage1>90)&(hsvImage1<=150));
        green1 = medfilt2(green1, [5,5]);
        binaryImage_G = bwareafilt(green1, count);
        labeledImage_G = bwlabel(binaryImage_G);
        % get information of the box [left, top, width, height]
        measurements_G = regionprops(labeledImage_G, 'BoundingBox');
        Pos_G = position(measurements_G, 'g', seg_width, box_size_target(1));
    
        % subtract blue colour with RGB
        blue = medfilt2(imsubtract(image(:,:,3), img_gray_target), [6,6]);
        blue1 = im2bw(blue,0.17); 
        binaryImage_B = bwareafilt(blue1, count);
        labeledImage_B = bwlabel(binaryImage_B);
        % get information of the box [left, top, width, height]
        measurements_B = regionprops(labeledImage_B, 'BoundingBox');
        Pos_B = position(measurements_B, 'b', seg_width, box_size_target(1));
    
        % subtract yellow colour with HSV
        yellow = ((hsvImage1>40)&(hsvImage1<=90));
        yellow1 = medfilt2(yellow, [5,5]);
        binaryImage_Y = bwareafilt(yellow1, count); 
        labeledImage_Y = bwlabel(binaryImage_Y);
        % get information of the box [left, top, width, height]
        measurements_Y = regionprops(labeledImage_Y, 'BoundingBox');
        Pos_Y = position(measurements_Y, 'y', seg_width, box_size_target(1));
    
        % subtract white colour using HSV
        whitePixels = s < 0.1 & v > 0;
        white1 = medfilt2(whitePixels, [8,8]);
        binaryImage_W = bwareafilt(white1, count);
        labeledImage_W = bwlabel(binaryImage_W);
        % get information of the box [left, top, width, height]
        measurements_W = regionprops(labeledImage_W, 'BoundingBox');
        %measurements_W = regionprops(binaryImage_W, 'BoundingBox');
        Pos_W = position(measurements_W, 'w', seg_width, box_size_target(1));
    
        % Combine all colours regions
        pos_combine = [Pos_R, Pos_G, Pos_B, Pos_Y, Pos_W];
        % pos_combine will keep information of Xmin, Ymin, colour code for
        % each colour regions
    
        % sorting the position
        pos_combine = sort(pos_combine);
    
        % Create matrix
        pos_x = 4;
        pre_result = [pos_combine(3)];    %extract the colour code for cumulative 3 order
        for item = 2:length(pos_combine)/3
            pre_result = [pre_result, pos_combine(pos_x + 2)];
            pos_x = pos_x + 3;
        end    
    
        % create 4x4 Cell matrix
        matrix_result = cell(rows,cols);
        g = 1;
        for r = 1:rows
            for c = 1:cols
                  matrix_result{r,c} = pre_result(g);  
                  if g < length(pre_result)
                    g = g + 1;
                end    
            end    
        end
    end
    
    % function to find circles  
    function [area, Centroid] = box(image)
        % filter the bw image
        n = 1;
        binaryImage = bwareafilt(im2bw(image,0.7), 1);
        binaryImage = medfilt2(binaryImage, [1 1]);
        % label the non black objects
        labeledImage = bwlabel(~binaryImage,8); % extract non black objects
        m = 1; 
        Centroid = [];
        % get region of all non black objects
        measurements = regionprops(labeledImage, 'BoundingBox', 'centroid');
        mea = measurements(1).BoundingBox(3);
        for i = 2:length(measurements)
            mea_1 = measurements(i).BoundingBox(3);
            if mea_1 > 10   % check size of the object to avoid taking noise
                if mea < mea_1
                    mea = mea_1;
                    n = i;
                else
                    Centroid = [Centroid; measurements(m).Centroid];
                    m = i;
                end
            end    
        end
        area = measurements(n).BoundingBox;
        Centroid = [Centroid; measurements(m).Centroid];
    end % end of box function

    % function to extract position of each segmented colour rectangles
    function result_seg = position(measure,colInit, limit, start)
        mea_arr = []
        for x = 1:length(measure)
            mr = measure(x).BoundingBox;
            % in case the rectangle is not properly separated, check how many
            % rectangle are joined together 
            diff = abs(mr(3) - mr(4));
            if ((mr(1) >= start) && (mr(3) >= limit - 20))
                    if mr(3) > mr(4)   % if they join horizontally
                        n = round(mr(3) / mr(4));
                        m1 = mr(1);
                        for z = 1:n
                            % return Xmin, Ymin and color code
                            mea_arr = [mea_arr, m1, mr(2), colInit];
                            m1 = m1 + mr(4);
                        end
                    else  % if they join vertically
                        n = round(mr(4) / mr(3));
                        m2 = mr(2);
                        for z = 1:n
                            % return Xmin, Ymin and color code
                            mea_arr = [mea_arr, mr(1), m2, colInit];
                            m2 = m2 + mr(3);
                        end
                    end    
            end   
        end
        result_seg = mea_arr;
    end % end of function result_seg
    
    % sort function
    function new_list = sort(list)
        n = (length(list)/3);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            m = 4;
            for i = 2:n
                % Swap elements in wrong order
                m_2 = list(m + 1);
                m_3 = list(m-2);
                if ((m_3 - m_2) > 20)
                    list = swap(list,m, m - 3);
                    nnew = i;
                elseif (m_2 <= m_3)
                    if (list(m) < list(m-3))
                        list = swap(list,m,m - 3);
                        nnew = i;
                    end
                elseif ((m_2 - m_3) < 20)
                    if (list(m) < list(m-3))
                        list = swap(list,m,m - 3);
                        nnew = i;
                    end    
                end
                m = m + 3;
            end
            n = nnew;
        end
        new_list = list;
    end % end of sort function   
    
    % function to swap object of sorting function 
    function x = swap(x,i,j)
    % Swap x(i) and x(j)
        val_1 = x(i);
        val_2 = x(i+1);
        val_3 = x(i+2);
        x(i) = x(j);
        x(i+1) = x(j+1);
        x(i+2) = x(j+2);
        x(j) = val_1;
        x(j+1) = val_2;
        x(j+2) = val_3;
    end % end of swap function
end

