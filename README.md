
# Point Cloud Alignment by Extremal Points
The main idea of this solution is to find the extreme points of the cloud:
- the furthest points in the cloud (YV and VB) 
- the furthest point (Base) from the YV - VB axis
- formation of transformation matrix for source cloud to target cloud.
- centroidal alignment
## Overview

This process involves aligning two point clouds based on three key points: the vertex (YV), the base point (YB), and an additional point farthest from the axis defined by YV and YB (referred to as "Base"). The alignment follows these steps:

1. **Align the vertices (YV)** of both models by adjusting their coordinates.
2. **Align the YB points** such that the YB of the source cloud lies on the YV-YB line of the target cloud.
3. **Align the Base points** by rotating the source cloud to match the plane formed by the YV, YB, and Base points of the target cloud.

## Initial Alignment (Preliminary/Trial)

The preliminary alignment between two point clouds is based on the three extremal points:

1. **YV (Vertex)**: The point furthest away from the base axis (largest distance from the other points).
2. **YB (Base)**: The point with the smaller distance to the "Base" point (which is farthest from the axis Y1-Y2).
3. **Base Point**: The point farthest from the axis defined by the YV and YB points.

In this preliminary step, we attempt to position the YV points of both clouds so they coincide, followed by aligning the YB point from the source onto the YV-YB line of the target, and ensuring the Base point aligns with the plane formed by the target YV, YB, and Base points.

## Alignment Procedure

1. **Align YV Points**
   - Calculate the translation needed to move the source YV to coincide with the target YV.

2. **Align YB Points**
   - Rotate the source cloud around the YV point to align the YB points on the same axis.

3. **Align Base Points**
   - Further rotate the source cloud around the YV-YB axis to align the Base point with the target Base point in the X-Z plane.

4. **Refine Alignment**
   - Fine-tune the transformation by ensuring the Source Base aligns with the plane defined by the Target YV, YB, and Base points.

## Steps for Optimization

- Improve efficiency in identifying extremal points (YV, YB, and Base).
- Refine rotations to ensure minimal deviations along the Z-axis.
- Perform tests with different point clouds to verify stability of the alignment procedure.

You can run this code form "start.bat" or from command line:
## python icp_script.py source_points.json target_points.json 0.001 50 0 1
### where:
- source_points.json  source cloude
- target_points.json  target cloude
- 0.001 tollerance
- 50 max itterations for centroidal alignment
- 0 - using algoritm 0 = by Extremal Points & centroidal alignment 
                     1 = by Extremal Points 
                     2 = by centroidal alignment
- 1 = visualise the result; 0 = without visualisation

## Features
Preliminary Alignment: The script identifies the key points (YV, YB, and Base) of both source and target point clouds to perform a rough alignment. The process includes translation and rotation to align key points.
ICP Fine-Tuning: Optionally, the script applies ICP to fine-tune the alignment by iteratively minimizing the distance between points in the source and target clouds.
Customizable Tolerance and Iteration Limits: The script supports custom settings for convergence tolerance and maximum iterations during the ICP phase.
## Dependencies
This script requires the following Python packages:

numpy
scipy
matplotlib
joblib
You can install the required packages using pip:

bash:
### pip install numpy scipy matplotlib joblib 

mailto:alekstar64@gmail.com
