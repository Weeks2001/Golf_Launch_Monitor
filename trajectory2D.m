%% Trajectory Calculation using Physics Equations
golfBallLocationTable = [647.74689, 308.21036;
 680.69934, 268.65640;
 750.92932, 205.32681;
 821.30310, 156.37511;
 864.62378, 92.646408];

% Clean up any potential NaN values
validIndices = all(~isnan(golfBallLocationTable), 2);
cleanData = golfBallLocationTable(validIndices, :);

% Flip the y-coordinates
frameHeight = 720; 
cleanData(:,2) = frameHeight - cleanData(:,2);

% Calculate time vector (assuming constant frame rate)
frameRate = 240; 
timeVector = (0:size(cleanData,1)-1)' / frameRate;

% Calculate velocity (pixels per second)
velocityX = diff(cleanData(:,1)) * frameRate;
velocityY = diff(cleanData(:,2)) * frameRate;
velocity = sqrt(velocityX.^2 + velocityY.^2);

% Calculate acceleration (pixels per second^2)
accelerationX = diff(velocityX) * frameRate;
accelerationY = diff(velocityY) * frameRate;
acceleration = sqrt(accelerationX.^2 + accelerationY.^2);

% Pixel to yards conversion factor 
ballDiameterInches = 1.68;
ballDiameterPixels = 125; 
inchesPerPixel = ballDiameterInches / ballDiameterPixels;
yardsPerPixel = inchesPerPixel / 36; % 36 inches = 1 yard

% Convert trajectory data to yards
xPositionYards = cleanData(:,1) * yardsPerPixel;
yPositionYards = cleanData(:,2) * yardsPerPixel;

vX_yards = velocityX * yardsPerPixel;
vY_yards = velocityY * yardsPerPixel;

% Ground level 
groundLevelYards = yPositionYards(1);

%% Physics-based calculations
% Get initial velocities (from the first frame)
initialVelocityX = velocityX(1) * yardsPerPixel; % yards per second
initialVelocityY = velocityY(1) * yardsPerPixel; % yards per second
initialSpeed = sqrt(initialVelocityX^2 + initialVelocityY^2); % yards per second
launchAngle = atan2d(initialVelocityY, initialVelocityX); % degrees

% Calculate acceleration due to gravity in yards/s^2
% Standard gravity is 9.81 m/s^2 = 32.2 ft/s^2 = 10.73 yards/s^2
gravity = 10.73; % yards per second squared

% Calculate time to reach maximum height (when vertical velocity becomes zero)
timeToMaxHeight = initialVelocityY / gravity;

% Calculate maximum height using equation: h = h₀ + v₀t - 0.5gt²
maxHeight = groundLevelYards + initialVelocityY * timeToMaxHeight - 0.5 * gravity * timeToMaxHeight^2;
maxHeightFromGround = maxHeight - groundLevelYards;

% Calculate time to return to ground level
% Use quadratic formula: h = h₀ + v₀t - 0.5gt²
% We want to find when h = h₀ (ground level)
% So: 0 = v₀t - 0.5gt²
% This simplifies to: t = 2v₀/g
timeToGround = 2 * initialVelocityY / gravity;

% Calculate final x position (total distance traveled)
finalXPositionYards = xPositionYards(1) + initialVelocityX * timeToGround;
totalDistanceYards = finalXPositionYards - xPositionYards(1);

% Generate physics-based trajectory points for plotting
tPhysics = linspace(0, timeToGround, 200)';
xPhysicsYards = xPositionYards(1) + initialVelocityX * tPhysics;
yPhysicsYards = groundLevelYards + initialVelocityY * tPhysics - 0.5 * gravity * tPhysics.^2;

% Create figures
% Position plot in yards
figure;
plot(xPositionYards, yPositionYards, 'bo-', 'LineWidth', 2);
hold on;
plot(xPhysicsYards, yPhysicsYards, 'r--', 'LineWidth', 1.5);
plot(finalXPositionYards, groundLevelYards, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
plot(xPositionYards(1) + initialVelocityX * timeToMaxHeight, maxHeight, 'mx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('X position (yards)');
ylabel('Y position (yards)');
title('Golf Ball Trajectory (Physics-Based)');
legend('Measured Data', 'Physics-Based Path', 'Predicted Landing Point', 'Maximum Height');
grid on;
axis equal;

% Position vs time
figure;
subplot(2,1,1);
plot(timeVector, xPositionYards, 'bo-');
hold on;
plot(tPhysics, xPhysicsYards, 'r--');
xlabel('Time (s)');
ylabel('X position (yards)');
title('X Position vs Time');
grid on;

subplot(2,1,2);
plot(timeVector, yPositionYards, 'bo-');
hold on;
plot(tPhysics, yPhysicsYards, 'r--');
plot([timeToGround timeToGround], [min(yPhysicsYards) max(yPhysicsYards)], 'g--');
plot([timeToMaxHeight timeToMaxHeight], [min(yPhysicsYards) max(yPhysicsYards)], 'm--');
xlabel('Time (s)');
ylabel('Y position (yards)');
title('Y Position vs Time');
grid on;
legend('Measured Data', 'Physics-Based Path', 'Ground Impact', 'Maximum Height');

% Convert velocities to mph (1 yard/s = 2.04545 mph)
initialVelocityMPH = initialSpeed * 2.04545;

% Display physics-based results
disp('Golf Ball Trajectory Analysis (Physics-Based):');
fprintf('Ball diameter used for scaling: %.2f pixels (%.2f inches)\n', ballDiameterPixels, ballDiameterInches);
fprintf('Conversion factor: %.6f yards/pixel\n', yardsPerPixel);
fprintf('Initial velocity: %.2f mph (%.2f yards/s)\n', initialVelocityMPH, initialSpeed);
fprintf('Launch angle: %.2f degrees\n', launchAngle);
fprintf('Estimated carry distance: %.2f yards\n', totalDistanceYards);
fprintf('Predicted landing position: %.2f yards\n', finalXPositionYards);
fprintf('Time in air: %.2f seconds\n', timeToGround);
fprintf('Maximum height: %.2f yards (reached at %.2f seconds)\n', maxHeightFromGround, timeToMaxHeight);

% 3D Trajectory Plot
zPhysicsYards = zeros(size(tPhysics));  
zMeasured = zeros(size(xPositionYards));  
vZ_yards = zeros(size(vX_yards)); 

figure('Position', [300 100 800 800]);
plot3(xPositionYards, zMeasured, yPositionYards, 'bo-', 'LineWidth', 2); hold on;
plot3(xPhysicsYards, zPhysicsYards, yPhysicsYards, 'r--', 'LineWidth', 1.5);
plot3(finalXPositionYards, 0, groundLevelYards, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
plot3(xPositionYards(1) + initialVelocityX * timeToMaxHeight, 0, maxHeight, 'mx', 'MarkerSize', 12, 'LineWidth', 2);

xlabel('X Position (yards)');
ylabel('Z (Lateral, yards)');
zlabel('Y Position (Height, yards)');
title('3D Golf Ball Trajectory');
grid on;
axis equal;
view(45, 25);

% Set axis limits
xlim([0 5]);
ylim([-2 2]);        
zlim([0 5]);        

% Midpoints for vector positions
x3_mid = x_mid;
y3_mid = zeros(size(x3_mid));  % Assuming no Z/lateral motion
z3_mid = y_mid;

figure('Position', [300 100 800 800]);
plot3(xPositionYards, zeros(size(xPositionYards)), yPositionYards, 'bo-', 'LineWidth', 2); hold on;