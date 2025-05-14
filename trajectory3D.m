% -----------------------------
% Step 1: Tracked Pixel Locations 
% -----------------------------
rawPixelCoords = [
    467.68359, 935.88629;
    467.12167, 935.75519;
    467.30258, 936.41687;
    467.30258, 936.41687;
    467.30258, 936.41687;
    640.58875, 573.28528;
    737.00000, 384.00000;
    805.00000, 229.00000;
    869.00000, 102.00000
];

% Remove stationary frames
pixelDiffs = sqrt(sum(diff(rawPixelCoords).^2, 2));
startIdx = find(pixelDiffs > 1, 1);
trackedPixels = rawPixelCoords(startIdx:end, :);

% Interpolate missing frames if there's a large jump
jumpThreshold = 100;
motionDiffs = sqrt(sum(diff(trackedPixels).^2, 2));
jumpIdx = find(motionDiffs > jumpThreshold);
if ~isempty(jumpIdx)
    i = jumpIdx(1);
    pt1 = trackedPixels(i, :);
    pt2 = trackedPixels(i+1, :);
    trackedPixels = [trackedPixels(1:i, :); ...
                     pt1 + (pt2-pt1)/3; ...
                     pt1 + 2*(pt2-pt1)/3; ...
                     trackedPixels(i+1:end, :)];
end

% -----------------------------
% Step 2: Pixel to Physical Conversion
% -----------------------------
ballDiameterInches = 1.68;
ballDiameterPixels = 90;
inchesPerPixel = ballDiameterInches / ballDiameterPixels;
yardsPerPixel = inchesPerPixel / 36;

% Flip Y-coordinates (origin at top-left)
frameHeight = 1080;
trackedPixels(:,2) = frameHeight - trackedPixels(:,2);

% Convert to yards
x_yards_measured = trackedPixels(:,1) * yardsPerPixel;
y_yards_measured = trackedPixels(:,2) * yardsPerPixel;

% -----------------------------
% Step 3: Compute v_x and v_y from pixels
% -----------------------------
fps = 240;
timeVector = (0:length(x_yards_measured)-1)' / fps;
v_x = diff(x_yards_measured) * fps;
v_y = diff(y_yards_measured) * fps;

% Use first frame's velocities as estimates
vx0 = v_x(1);
vy0 = v_y(1);

% -----------------------------
% Step 4: Solve for v_z using known total speed
% -----------------------------
initialSpeedMPH = 152.4;
v0 = initialSpeedMPH / 2.04545;  % total speed in yards/sec
vz0 = sqrt(v0^2 - vx0^2 - vy0^2);  % derived forward speed

% -----------------------------
% Step 4b: Calculate Launch Direction (Azimuth Angle)
% -----------------------------
launchAzimuthDeg = atan2d(vx0, vz0);  % angle relative to Z-axis

fprintf('Launch direction (azimuth): %.2f degrees\n', launchAzimuthDeg);

% -----------------------------
% Step 5: Simulate 3D Trajectory
% -----------------------------
g = 10.73;
timeToGround = 2 * vy0 / g;
t = linspace(0, timeToGround, 300)';

x_traj = vx0 * t;
z_traj = vz0 * t;
y_traj = vy0 * t - 0.5 * g * t.^2;

[maxY, maxIdx] = max(y_traj);

% -----------------------------
% Step 6: Plot 3D Trajectory
% -----------------------------
figure('Position', [200 200 900 600]);
plot3(x_traj, z_traj, y_traj, 'r-', 'LineWidth', 2); hold on;
plot3(x_traj(maxIdx), z_traj(maxIdx), maxY, 'bx', 'MarkerSize', 12, 'LineWidth', 2);
plot3(x_traj(end), z_traj(end), 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2);

xlabel('Lateral X (yards)');
ylabel('Forward Z (yards)');
zlabel('Height Y (yards)');
title(sprintf('3D Golf Ball Trajectory (%.1f mph)', initialSpeedMPH));
legend('Trajectory', 'Max Height', 'Landing Point');
grid on;
xlim([-25, 25]);
ylim([0, 300]);
zlim([0, 20]);
view(45, 25);

% -----------------------------
% Step 7: Display Output
% -----------------------------
fprintf('Estimated v_x from pixels: %.2f yards/sec\n', vx0);
fprintf('Estimated v_y from pixels: %.2f yards/sec\n', vy0);
fprintf('Derived v_z: %.2f yards/sec\n', vz0);
fprintf('Initial speed: %.2f mph (%.2f yards/sec)\n', initialSpeedMPH, v0);
fprintf('Total forward distance: %.2f yards\n', z_traj(end));
fprintf('Maximum height: %.2f yards\n', maxY);
fprintf('Time in air: %.2f seconds\n', timeToGround);
