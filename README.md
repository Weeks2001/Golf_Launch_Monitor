Project Title
Predicting Golf Ball Trajectory Using Short Burst Imagery of Swing

Authors
Luke Bagan
Colby J Weeks

Course
CSCE 867: Computer Vision
University of South Carolina, Columbia
May 6, 2025

Overview
This project presents a low-cost, accessible golf swing analysis system that estimates golf ball trajectory using high-frame-rate consumer cameras and standard laptops. By combining core computer vision techniques with physics-based modeling, the system predicts real-world golf ball flight characteristics like launch angle, velocity, maximum height, and distance traveled. The goal is to offer an affordable alternative to commercial tracking systems, giving golfers valuable insights into their swing mechanics without specialized hardware.

Key Features
Camera calibration using checkerboard patterns for 2D-3D correspondence
Ball detection post-impact using lightweight tracking algorithms
Frame-by-frame motion tracking and pixel-to-world conversion
Trajectory prediction using classical kinematic equations
Comparison against ground-truth data for accuracy assessment
2D proof-of-concept and 3D trajectory modeling integration

Methodology
Ball Detection
High-frame-rate video footage is analyzed to detect the golf ball after club contact.
A simple feature-based tracking algorithm is employed to follow the ball's motion frame-by-frame.
Pixel coordinates are converted to real-world distances using the known size of a golf ball.
Trajectory Calculation
Velocity, acceleration, and launch angle are calculated from frame-to-frame position changes.
Standard kinematic equations are used to estimate the ball's flight path.
Parameters like maximum height, total distance, and flight time are predicted.

Modeling & Test Setup
2D Testing: GoPro camera at 240 FPS (rolling shutter) for short-range analysis.
3D Testing: Rapsodo launch monitor with a global shutter for higher fidelity data.
Launch monitor speed data was combined with visual tracking for full 3D trajectory modeling.

Results
Predicted distance within 15.8% error of actual distance.
Launch direction estimated with 12% error.
Maximum height estimation had larger error (~69%) due to camera angle distortion.
Demonstrated feasibility of vision-based golf ball tracking using consumer-grade equipment.

Challenges
Motion blur from rolling shutter cameras
Occlusion of the ball by the club during early frames
Depth estimation without stereo vision
Pixel-to-world scaling inaccuracies due to camera geometry

Future Work
Integrating stereo vision for true 3D positional tracking
Improving detection algorithms (e.g., template matching, AI-based models)
Implementing Kalman filtering for smoother tracking
Adding spin analysis and environmental factors (drag, wind, altitude)

Potential integration of radar or LiDAR for enhanced accuracy
