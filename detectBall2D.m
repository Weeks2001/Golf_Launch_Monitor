clear all
close all

folderPath = 'C:\Users\lbaga\Desktop\CSCE 867\Final Project\highSpeed';
folder = dir(folderPath);
folder = folder(3:end);
video = VideoWriter('highSpeedDetection.mp4');
open(video)
for i = 1:length(folder)
    I = imread(string(folderPath) + '\' + string(folder(i).name));

    Mask = rgb2gray(I)>220;
    Mask = imopen(Mask, strel('disk',10));
    Mask = imclose(Mask, strel('disk',12));

    [regions,cc] = detectMSERFeatures(Mask, 'RegionAreaRange',[1000,size(I,1)/4*size(I,2)/4], 'ThresholdDelta',4);
    stats = regionprops('table',cc,'Eccentricity');
    eccentricityIdx = stats.Eccentricity > 0.6;
    nonCircularRegions = regions(eccentricityIdx);
    [nonCircFeat, valPoints] = extractFeatures(I(:,:,3), nonCircularRegions);

    figure
    imshow(I); hold on;
    plot(nonCircularRegions,'showPixelList',false,'showEllipses',true);

    hasGolfBallLeftFrame = 0;

    if length(nonCircularRegions) == 1
        if i ~=1
            if numFeatPrev == 2
                hasGolfBallLeftFrame = 1;
                leftFrameidx = i;
            end
        end
        golfBallFeat = nonCircFeat;
        golfBallPoints = valPoints;
        golfBallLocation = valPoints.Location;
    elseif length(nonCircularRegions) > 1
        idxPairs = matchFeatures(golfBallFeat,nonCircFeat);
        golfBallFeatMatch = valPoints(idxPairs(:,1));
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
            golfBallPoints = valPoints;
            golfBallLocation = golfBallFeatMatch.Location;
        end

    end

    numFeatPrev = length(nonCircularRegions);

    golfBallLocationTable(i,:) = golfBallLocation;

    hold on; plot(golfBallLocation(1), golfBallLocation(2),'ro')

    F = getframe(gcf);
    [X, Map] = frame2im(F);
    writeVideo(video,X);

end

golfBallLocationTable(leftFrameidx:end,:) = [];
close(video)