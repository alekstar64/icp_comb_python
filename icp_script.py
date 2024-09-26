# -*- coding: utf-8 -*-
from scipy.spatial import KDTree
import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import pdist, squareform
from joblib import Parallel, delayed
#from scipy.spatial import KDTree
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#from trimesh import transform_points

def transformed_points(points,transformation):
    points_homogeneous = np.hstack((points, np.ones((points.shape[0], 1))))
    transformation_matrix = np.array(transformation)
    points_transformed_homogeneous = points_homogeneous @ transformation_matrix.T
    points_transformed = points_transformed_homogeneous[:, :3]
    return points_transformed
def distance_points(point1, point2):
    return np.linalg.norm(point2 - point1)
def thin_out_points(points, N):
    return points[::N]
def find_extreme_points(points):
    """
    Find the two farthest points (Y1 and Y2) in the cloud and the point farthest from the line Y1-Y2 (Base).

    Parameters:
    points (numpy array): Nx3 array of 3D points.

    Returns:
    tuple: (YV, YB, Base) where YV and YB are the two farthest points and Base is farthest from YV-YB line.
    """
    if len(points) > 1000:
        points = points[::5]
    distances = pdist(points)
    distance_matrix = squareform(distances)
    i, j = np.unravel_index(np.argmax(distance_matrix), distance_matrix.shape)
    Y1, Y2 = points[i], points[j]

    def point_to_line_distance(point, line_point1, line_point2):
        """
        Calculate distance from a point to a line (defined by two points).
        """
        line_vec = line_point2 - line_point1
        point_vec = point - line_point1
        line_len = np.linalg.norm(line_vec)
        return np.linalg.norm(np.cross(line_vec, point_vec)) / line_len if line_len > 0 else np.linalg.norm(point_vec)

    max_base_dist = 0
    Base = points[0]
    for point in points:
        dist = point_to_line_distance(point, Y1, Y2)
        if dist > max_base_dist:
            max_base_dist = dist
            Base = point

    dist_Y1_to_Base = np.linalg.norm(Y1 - Base)
    dist_Y2_to_Base = np.linalg.norm(Y2 - Base)

    if dist_Y1_to_Base > dist_Y2_to_Base:
        YV, YB = Y1, Y2
    else:
        YV, YB = Y2, Y1

    return YV, YB, Base
def angle_between_vectors(v1, v2):
    # Normalize the vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    # Compute the dot product
    dot_product = np.dot(v1, v2)
    # Clamp the value to avoid numerical issues
    dot_product = np.clip(dot_product, -1.0, 1.0)
    # Compute the angle
    angle = np.arccos(dot_product)
    return angle

def rotation_matrix(axis, theta):
    # Normalize the axis
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)  
    # Compute quaternion components
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    a, b, c, d = [x.item() for x in [a, b, c, d]]
    # Construct the rotation matrix
    return np.array([
        [1 - 2 * (c**2 + d**2), 2 * (b*c - a*d), 2 * (b*d + a*c), 0],
        [2 * (b*c + a*d), 1 - 2 * (b**2 + d**2), 2 * (c*d - a*b), 0],
        [2 * (b*d - a*c), 2 * (c*d + a*b), 1 - 2 * (b**2 + c**2), 0],
        [0, 0, 0, 1]
    ], dtype=float)

# Helper function to find the perpendicular point of 'Point' on the line defined by YV-YB
def perpendicular_point_on_line(YV, YB, Point):
    line_vector = YB - YV
    line_vector_norm = line_vector / np.linalg.norm(line_vector)
    point_vector = Point - YV
    perpendicular_point = YV + np.dot(point_vector, line_vector_norm) * line_vector_norm
    return perpendicular_point
def align_point_clouds(target_YV, target_YB, target_Base, source_YV, source_YB, source_Base, trans_step):
    """
    Align the source point cloud to the target point cloud by first translating and then rotating the source points.
    
    Parameters:
    target_YV, target_YB, target_Base (numpy arrays): Target points.
    source_YV, source_YB, source_Base (numpy arrays): Source points.

    Returns:
    list: The final transformation matrix (4x4).
    """
    
    # Step 1: Translation to align YB points
    source_orig = np.array([source_YV, source_YB, source_Base])
    translation_to_origin = target_YB - source_YB
    translation_matrix_to_origin = np.eye(4)
    translation_matrix_to_origin[:3, 3] = translation_to_origin
    # Translate source points to new origin
    translated_source_YV = source_YV + translation_to_origin
    translated_source_YB = source_YB + translation_to_origin
    translated_source_Base = source_Base + translation_to_origin    
    # Step 2: Rotation to align YV points
    # Compute vectors for rotation
    vector_target = target_YV - target_YB
    vector_source = translated_source_YV - translated_source_YB  
    # Compute the rotation angle and axis
    angle_1 = angle_between_vectors(vector_target, vector_source)
    axis = np.cross(vector_target, vector_source)   
    axis = np.asarray(axis).flatten()   
    # Default axis if zero length
    if np.linalg.norm(axis) > 1e-6:
        axis /= np.linalg.norm(axis)
    else:
        axis = np.array([0, 0, 1])
    # Compute the rotation matrix around the origin
    rotation_matrix_xz = rotation_matrix(axis, angle_1)
    #print(f"Rotation matrix: {rotation_matrix_xz}")
    # Translate the source points so that target_YB becomes the origin
    source_trans_xz = transformed_points(source_orig, rotation_matrix_xz)
    translation_to_xz = target_YB - source_trans_xz[1]
    translation_matrix_to_xz = np.eye(4)
    translation_matrix_to_xz[:3, 3] = translation_to_xz
  
    # Translate source points to new origin
    source_trans_xz[0] = source_trans_xz[0] + translation_to_xz
    source_trans_xz[1] = source_trans_xz[1] + translation_to_xz
    source_trans_xz[2] = source_trans_xz[2] + translation_to_xz
   
    # Step 3: Rotate in the X-Y plane to align the YB points
    axis_yv_yb = target_YB - target_YV
    axis_yv_yb /= np.linalg.norm(axis_yv_yb)
    source_base_projection = perpendicular_point_on_line(target_YV, target_YB, source_trans_xz[2])
    target_base_projection = perpendicular_point_on_line(target_YV, target_YB, target_Base)
    source_base_projection = source_trans_xz[2] - source_base_projection 
    target_base_projection = target_Base - target_base_projection
    angle = -angle_between_vectors(source_base_projection, target_base_projection)

    # # Step 3: Rotate in the X-Y plane to align the YB points
    print("angle_XY" , angle)
    rotation_matrix_xy = rotation_matrix(axis_yv_yb, angle)    
    source_trans_x1 = transformed_points(source_orig,rotation_matrix_xy @ rotation_matrix_xz)
    source_trans_xy = transformed_points(source_orig,rotation_matrix_xy @ rotation_matrix_xz @ translation_matrix_to_origin)
    translation_matrix_to_xy = np.eye(4)
    translation_matrix_to_xy[:3, 3] = target_YB - source_trans_xy[1]
    final_transformation_matrix = (translation_matrix_to_xy @ rotation_matrix_xy @ rotation_matrix_xz @ translation_matrix_to_origin).tolist()
    return final_transformation_matrix

def find_key_points(points, tolerance=1e-5):
    # Step 1: Identify Y1 (Apex) as the point with the maximum Y value
    y1 = points[np.argmax(points[:, 1])]

    # Step 2: Identify Y2 as the base point farthest from Y1 in the X-Z plane
    base_points = points[np.abs(points[:, 1]) < tolerance]
    if base_points.shape[0] == 0:
        raise ValueError("No base points found with Y close to 0")
    max_distance = 0
    y2 = base_points[0]
    for point in base_points:
        distance = np.linalg.norm(point - y1)
        if distance > max_distance:
            max_distance = distance
            y2 = point
    # Step 3: Find the opposite point to Y2 in the X-Z plane
    max_distance = 0
    opposite_point = y2
    for point in base_points:
        distance = np.linalg.norm(point[[0, 2]] - y2[[0, 2]])  # Only compare X-Z coordinates
        if distance > max_distance:
            max_distance = distance
            opposite_point = point
    return y1, y2, opposite_point


def icp(source_points, target_points,max_iterations=50, tolerance=0.1):    
    source_points = np.array(source_points)
    target_points = np.array(target_points)    
    transformation = np.eye(4)
    
    # Refine the transformation with ICP if required
    for i in range(max_iterations):
        tree = KDTree(target_points)
        distances, indices = tree.query(source_points)
        closest_points = target_points[indices]

        source_centroid = np.mean(source_points, axis=0)
        target_centroid = np.mean(closest_points, axis=0)

        source_centered = source_points - source_centroid
        target_centered = closest_points - target_centroid

        H = np.dot(source_centered.T, target_centered)
        U, _, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)

        if np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = np.dot(Vt.T, U.T)

        t = target_centroid - np.dot(R, source_centroid)

        new_transformation = np.eye(4)
        new_transformation[:3, :3] = R
        new_transformation[:3, 3] = t

        transformation = np.dot(new_transformation, transformation)

        source_points = np.dot(source_points, R.T) + t

        if np.linalg.norm(t) < tolerance:
            #print(f"ICP converged after {i+1} iterations.")
            break
    print(f"ICP converged after {i+1} iterations.")
    return transformation.tolist()

if __name__ == "__main__":

    # axis = np.array([0, 0, 1])
    # theta = np.pi / 2  

    # rotation_mat = rotation_matrix(axis, theta)
    # print(f"Rotation matrix for 90 degrees around Z-axis:\n{rotation_mat}")


    # source_points = np.array([[1, 0, 0], [0, 1, 0]])
    # transformed_points = np.dot(rotation_mat, np.vstack((source_points.T, np.ones(source_points.shape[0]))))

    # print(f"Transformed points:\n{transformed_points[:3, :].T}")
    # exit


    source_file = sys.argv[1]
    target_file = sys.argv[2]    
    if len(sys.argv) > 3:
        tolerance = float(sys.argv[3])
    else:
        tolerance = 0.01
    if len(sys.argv) > 4:
        max_iterations = int(sys.argv[4])
    else:
        max_iterations = 20
    if len(sys.argv) > 5:
        scema = int(sys.argv[5])
    else:
        scema = 0
    if len(sys.argv) > 6:
        output = int(sys.argv[6])
    else:
        output = 0 # not show result
        
    # scema = 0 - find_extreme_points & align_point_clouds & icp
    # scema = 1 - only find_extreme_points & align_point_clouds        
    # scema = 2 - only icp
    with open(source_file, 'r') as f:
        source_points = np.array(json.load(f), dtype=np.float32)

    with open(target_file, 'r') as f:
        target_points = np.array(json.load(f), dtype=np.float32)
    target_YV = target_YB = target_Base = source_YV = source_YB = source_Base = 0.0
    if scema < 2:
        target_YV, target_YB, target_Base = find_extreme_points(target_points)
        source_YV, source_YB, source_Base = find_extreme_points(source_points)
    #print(f"distance target_YB - target_Base: {distance_points(target_YB, target_Base)},\ndistance source_YB - source_Base: {distance_points(source_YB, source_Base)}")
    #source_3d = np.array([source_YV, source_YB, source_Base])
    # # # # # # target_points = np.array([
    # # # # # # # # #     #[ 235.,  130., -264.],[ 240.,  100., -269.], [ 230.,  105., -259.]])
    # # # # # #      [ 1., 1., 1.],[11.,11.,1.],[ 11., 1., 1.]])
    # # # # # # # # #     #[ 0., 0., 0.],[10.,  2., 0.]])
    # # # # # # source_points = np.array([
    # # # # # # # # #     [0.,0.,10.],[0.,0.,0.],[-10.,0.,0.]])
    # # # # # # # # #     #[221.92522,120.00881,-245.86998],[236.20018,105.093414,-268.7557],[231.62068,98.24139,-256.22263]
    # # # # # #      [0.,0.,10.],[7.071,7.071,0.],[0.,0.,0.]])
    # # # # # # # # #     #[0.,0.,0.],[0.,10.,0.]])
    # # # # # # # # #target_YV = np.array([ 235.0,  130.0, -264.0])
    # # # # # # # # # #target_YB = np.array([240.0,  100.0, -269.0])
    # # # # # # # # # #target_Base = np.array([230.,  105., -259.])
    # # # # # # target_YV = np.array(target_points[0])
    # # # # # # target_YB = np.array(target_points[1])
    # # # # # # target_Base = np.array(target_points[2])

    # # # # # # source_YV = np.array(source_points[0])
    # # # # # # source_YB = np.array(source_points[1])
    # # # # # # source_Base = np.array(source_points[2])
    # # # # # # print(f"Target YV: {target_YV}, YB: {target_YB}, Base: {target_Base}")
    # # # # # # print(f"Source YV: {source_YV}, YB: {source_YB}, Base: {source_Base}")

    # # # #source_YV = np.array([221.92522,120.00881,-245.86998])
    # # # #source_YB = np.array([236.20018,105.093414,-268.7557])
    # # # #source_Base = np.array([231.62068,98.24139,-256.22263])
    # # # #exit
    # # # # transformation = align_point_clouds(target_YV, target_YB, target_Base, source_YV, source_YB, source_Base, 1)
    # # # # print(f"transformation :", transformation)
    # # # # # Convert source points to homogeneous coordinates by adding a column of 1s
    # # # # source_points_homogeneous = np.hstack((source_points, np.ones((source_points.shape[0], 1))))
    # # # # # Perform the transformation (4x4 matrix multiplication)
    # # # # transformation_matrix = np.array(transformation)
    # # # # source_points_transformed_homogeneous = source_points_homogeneous @ transformation_matrix.T
    # # # # source_points_transformed = source_points_transformed_homogeneous[:, :3]
    # # # # print(f"Transformed source translatited:\n{source_points_transformed}")

    #transformation = align_point_clouds(target_YV, target_YB, target_Base, source_YV, source_YB, source_Base, 1)
    #print(f"transformation :", transformation)
    
    #print(f"target_YV - source_YV: {target_YV - source_points_transformed[0]}, \ntarget_YB - transformed_source_YB: {target_YB - source_points_transformed[1]}")
    #print(f"target_Base - source_Base: {target_Base - source_points_transformed[2]}")
    #print(f"distance target_YV - transformed_source_YV: {distance_points(target_YV, source_points_transformed[0])},\ndistance target_YB - transformed_source_YB: {distance_points(target_YB, source_points_transformed[1])}")
    #print(f"distance target_Base - transformed_source_Base: {distance_points(target_Base, source_points_transformed[2])}")

    #print(f"source_points_transformed: {source_points_transformed}")

    #print(f"target_points: {target_points}")
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(target_points[:, 0], target_points[:, 1], target_points[:, 2], c='r', label='target_points')  
    # #ax.scatter(source_points_transformed[:, 0], source_points_transformed[:, 1], source_points_transformed[:, 2], c='b', label='source_points_transformed')
    # ax.scatter(source_points[:, 0], source_points[:, 1], source_points[:, 2], c='b', label='source_points')

    
    # ax.set_xlabel('X Axis')
    # ax.set_ylabel('Y Axis')
    # ax.set_zlabel('Z Axis')


    # ax.legend()


    #plt.show()
    #exit
    # Print transformed points
#    print(f"Transformed source points:\n{source_points_transformed}")
    transformation = np.eye(4)         
    full_transformation = np.eye(4)         
    icp_transformation = np.eye(4)     
    source_points_transformed = []
    if scema < 2:
        transformation = align_point_clouds(target_YV, target_YB, target_Base, source_YV, source_YB, source_Base, 1)        
        source_points_transformed = transformed_points(source_points,transformation)
    if scema == 0 :
        icp_transformation = icp(source_points_transformed, target_points,max_iterations, tolerance)        
    if scema == 2:
        icp_transformation = icp(source_points, target_points,max_iterations, tolerance)        
        full_transformation = icp_transformation
    if scema == 0:
        full_transformation = np.dot(icp_transformation, transformation).tolist()
    if scema == 1:
        full_transformation = transformation
    source_points_transformed = transformed_points(source_points,full_transformation)        
    if output != 0:
        fig = plt.figure()        
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(target_points[:, 0], target_points[:, 1], target_points[:, 2], c='r', label='target_points')
        ax.scatter(source_points_transformed[:, 0], source_points_transformed[:, 1], source_points_transformed[:, 2], c='b', label='source_points_transformed')
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        ax.legend()    
        plt.show()
    transformation_file = os.path.join(os.path.dirname(source_file), 'transformation.json')    
#    print("full_transformation",full_transformation)
    with open(transformation_file, 'w') as f:
        json.dump(full_transformation, f)
