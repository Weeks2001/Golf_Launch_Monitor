clear all
close all

folderPath = 'C:\Users\lbaga\Desktop\CSCE 867\Final Project\252yds';
folder = dir(folderPath);
folder = folder(3:end);
folder = natsortfiles(folder);
firstOcclusion = 0;
video = VideoWriter('252Detection.avi');
open(video)
for i = 1:length(folder)

    I = imread(string(folderPath) + '\' + string(folder(i).name));
    I = I(1:1250, 250:1250,:);
    I = rgb2gray(I);
    I = MSRetinex2(im2double(I), [5, 35, 150], [.1 .1], 8);
    I = uint8(I);
    [rows, columns] = size(I);

    background = imopen(I, strel('disk', 10));
    J = imsubtract(I, background);
    J = histeq(J);
    J = imgaussfilt(J, 2, "FilterSize", [7 7]);
    counts = histcounts(J,255);
    T = otsuthresh(counts);
    T = round(T*255);

    Mask = J<T;
    Mask = imopen(Mask, strel('disk',10));
    Mask = imclose(Mask, strel('disk',12));

    [regions,cc] = detectMSERFeatures(Mask, 'RegionAreaRange',[1000,round(size(I,1)/4*size(I,2)/4)], 'ThresholdDelta',4);
    stats = regionprops('table',cc,'Eccentricity');
    eccentricityIdx = stats.Eccentricity < 0.65;
    circularRegions = regions(eccentricityIdx);

    bottomHalf = circularRegions.Location(:,2) > size(I,1)/2;
    circularRegions = circularRegions(bottomHalf);

    [circFeat, valPoints] = extractFeatures(I, circularRegions);

    figure
    imshow(I); hold on;
    plot(circularRegions,'showPixelList',false,'showEllipses',true);

    if ~isempty(circularRegions)

        if length(circularRegions) == 1
            golfBallFeat = circFeat;
            golfBallPoints = valPoints;
            golfBallLocation = valPoints.Location;
        elseif length(circularRegions) > 1
            idxPairs = matchFeatures(golfBallFeat,circFeat);
            golfBallFeatMatch = valPoints(idxPairs(:,2));
            if isempty(golfBallFeatMatch)
                motion = valPoints.Location - golfBallLocation;
                prevMotion = golfBallLocation - golfBallLocationTable(i-2,:);
                prevMotion = repmat(prevMotion,[size(motion,1),1]);
                cosineSimilarity = [];
     
                for j = 1:length(prevMotion)
                    cosineSimilarity(j) = dot(motion(j,:),prevMotion(j,:))./(sqrt(sum(motion(j,:).^2)).*sqrt(sum(prevMotion(j,:).^2)));
                end
    
                [~,maxIdx] = max(cosineSimilarity);
    
                golfBallLocation = valPoints(maxIdx).Location;
            else
                golfBallLocation = golfBallFeatMatch.Location;
            end
        end
    else
        if firstOcclusion == 0 || firstOcclusion == 1
            firstOcclusion = firstOcclusion +1;
            occluddedIdx = i;
        else

            [regions,cc] = detectMSERFeatures(J, 'RegionAreaRange',[2000,round(size(I,1)/4*size(I,2)/4)], 'ThresholdDelta',4);
            stats = regionprops('table',cc,'Eccentricity');
            eccentricityIdx = stats.Eccentricity < 0.65;
            circularRegions = regions(eccentricityIdx);
            [circFeat, valPoints] = extractFeatures(I, circularRegions);

            if isempty(circFeat)

                numTargets = 10;
                level = 1;
                threshold = double(.6);
                tempMean = mean2(originalTemplate);
                template = originalTemplate - tempMean;
                templateEnergy = sqrt(sum(template(:).^2));
                templateRotated = imrotate(template, 180);
                [tempR, tempC] = size(templateRotated);
                r_mod = 2^nextpow2(tempR + rows);
                c_mod = 2^nextpow2(tempC+columns);
                template_p = [templateRotated zeros(tempR,c_mod-tempC)];
                template_p = [template_p; zeros(r_mod-tempR, c_mod)];

                template_fft = fft2(template_p);
                target_size = repmat(size(template), [numTargets,1]);

                gain = 2^(level-1);
                Im_p = zeros(r_mod, c_mod, 'double');
                c_ones = ones(tempR, tempC, 'double');

                hFindMax = vision.LocalMaximaFinder('Threshold', double(-1), 'MaximumNumLocalMaxima', numTargets, 'NeighborhoodSize', floor(size(template)/2)*2-1);

                Im_p(1:rows, 1:columns) = (I-tempMean);
                I2_fft = fft2(Im_p);
                corr_freq = I2_fft .* template_fft;
                corr_out_f = ifft2(corr_freq);
                corr_out_f = corr_out_f(tempR:rows, tempC:columns);

                IUT_energy = (I-tempMean).^2;
                IUT = conv2(IUT_energy, c_ones, 'valid');
                IUT = sqrt(IUT);

                norm_corr_f = (corr_out_f)./(IUT.*templateEnergy);

                if size(norm_corr_f,1) < size(template,1) || size(norm_corr_f,2) < size(template,2)
                    continue
                end
                xyLocation = step(hFindMax, norm_corr_f);
                linear_index = sub2ind([rows-tempR, columns-tempC]+1, xyLocation(:,2), xyLocation(:,1));
                norm_corr_f_linear = norm_corr_f(:);
                norm_corr_val = norm_corr_f_linear(linear_index);
                detect = (norm_corr_val > threshold);

                percent_corr_val_dropoff = (norm_corr_val(1:end-1)-norm_corr_val(2:end))./norm_corr_val(1:end-1);
                largeDrops = percent_corr_val_dropoff > .05;

                if sum(largeDrops) > 0
                    detect(logical([largeDrops;0])) = 1;
                end
                target_roi = zeros(length(detect), 4);
                ul_corner = (gain.*(xyLocation(detect,:) -1))+1;

                target_roi(detect, :) = [ul_corner, fliplr(target_size(detect, :))];
                target_roi(target_roi(:,4) == 0, :) = [];

                searchArea = [golfBallLocation(1) golfBallLocation(2)-325 300 300];

                overlap = bboxOverlapRatio(target_roi, searchArea, 'Min');
                templateMatched = overlap>=.7;
                targetCentroids = [target_roi(:,1)+round((target_roi(:,3)/2)) target_roi(:,2)+round((target_roi(:,3)/2))];
                golfBallLocation = targetCentroids(templateMatched,:);
                if isempty(golfBallLocation)
                    golfBallOutOfFrameIdx = i;
                    break
                end


            else
            
                euclDist = sqrt(sum((circFeat-golfBallFeat).^2,2));
                [~,closestMatch] = min(euclDist);
                golfBallFeatMatch = valPoints(closestMatch);
                golfBallLocation = golfBallFeatMatch.Location;
                originalTemplate = I(golfBallLocation(2)-90:golfBallLocation(2)+90,golfBallLocation(1)-90:golfBallLocation(1)+90);
                hold on; plot(golfBallFeatMatch)

            end

        end

    end

    numFeatPrev = length(circularRegions);

    golfBallLocationTable(i,:) = golfBallLocation;

    hold on; plot(golfBallLocation(1), golfBallLocation(2),'ro')

    F = getframe(gcf);
    [X, Map] = frame2im(F);
    writeVideo(video,X);

end

golfBallLocationTable(golfBallOutOfFrameIdx:end,:) = [];
close(video)