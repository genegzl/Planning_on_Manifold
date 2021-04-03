#ifndef  MOTION_MANIFOLD_HPP
#define  MOTION_MANIFOLD_HPP

#include <Eigen/Dense>
// #include <ikd-Tree/ikd_Tree.h>
#include <ikd-Forest/ikd_Forest.h>
#include <ros/ros.h>
#include "std_msgs/Bool.h"
#include "std_msgs/Float32MultiArray.h"
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>

#define print_line printf("%s -- %d\r\n", __FILE__, __LINE__);

// using namespace std;
using namespace Eigen;


#define TRAJ_DATA_ROW                       1759
#define TRAJ_DATA_COL                       18
#define TRAJ_REPEAT_TICK                    0
#define ESTI_PLANE_PONIT_NUM                5
#define SURFACE_ESTI_PONIT_NUM_MIN          10
#define SURFACE_ESTI_PONIT_NUM_MAX          200
#define MIN_MAP_SIZE                        2.5e4
#define SURFACE_FITTING_THRESHOLE           100
#define SURFACE_SEARCH_RANGE_FORWARD        1.5
#define SURFACE_SEARCH_RANGE_BACK           1.5
#define SURFACE_SEARCH_RANGE_LEFT           1.5
#define SURFACE_SEARCH_RANGE_RIGHT          1.5
#define SURFACE_SEARCH_RANGE_UP             1.5
#define SURFACE_SEARCH_RANGE_DOWN           1.5
#define LIDAR_MOUNT_HIGHT                   0.3

extern std_msgs::Bool traj_request;

typedef float traj_data_t[TRAJ_DATA_ROW][TRAJ_DATA_COL];

enum traj_mode_t  
{
    HOVER           = 10,
    TRAJ_REPEAT     = 11,
    TRAJ_STAY       = 12
};

struct coordinate_t
{
    Vector2f    pos;
    Vector2f    vel;
    Vector2f    acc;
    float       yaw;
    float       yaw_rate;
    float       vb;
    float       ab;
};

struct traj_t
{
    bool                                               init_en;
    traj_mode_t                                        mode;
    uint32_t                                           wait_cnt;
    uint32_t                                           read_cnt;

    coordinate_t                                       coordinate;

    Vector3f                                           pos;                     // local position,      unit: m
    Vector3f                                           vel;                     // local velocity,      unit: m/s
    Vector3f                                           acc;                     // local acceleration,  unit: m/s^2
    Vector3f                                           omega;                   // body angular rate,   unit: rad/s
    float                                              yaw;                     // Euler yaw,           unit: rad
    Quaternionf                                        q;                       // quaternion,          order: wxyz
    Matrix3f                                           R;
    float                                              vb;
    float                                              ab;

    Quaternionf                                        q_init;                  // initial attitude   
    Vector3f                                           pos_init;
    float                                              yaw_init; 
 };

 struct plane_t
 {
     Vector3f                                   normal_vector;
     Vector3f                                   centroid;
     Vector3f                                   nearest_point;
     PointVector                                neighbor_ponits;
     MatrixXf                                   lsm_A;
     MatrixXf                                   lsm_b;
 };

 struct quadradic_surface_t
 {
     Matrix<float, 6, 1>                            para;
     PointVector                                    neighbor_ponits;
     uint32_t                                       neighbor_point_num;
     vector<float>                                  Distance2Centrois;
     Vector3f                                       search_center;
     BoxPointType                                   search_box;
     Vector3f                                       normal_vector;
     float                                          esti_err;
     float                                          p_F_x;
     float                                          p_F_y;
     float                                          p_F_xx;
     float                                          p_F_yy;
     float                                          p_F_xy;
     Vector3f                                       x_proj;
     Vector3f                                       y_proj;
     MatrixXf                                       lsm_A;
     MatrixXf                                       lsm_b;
 };

 struct ekf_data_t
 {
     Matrix3f       R;
     Vector3f       pos;
     Vector3f       vel;
 };

 struct traj_manifold_t
 {
    traj_t                          traj;
    plane_t                         tangent_plane;
 };

class MANIFOLD
{
private:
    bool                            init_en;
    bool                            surfaceValid;
    
    bool plane_estimate(plane_t &plane, float threshold, int point_num);
    void get_centroid(Vector3f &centroid, PointVector &point, int point_num);
    void point_prejoct_on_plane(Vector3f &p_out, Vector3f &p_in, plane_t &plane);
    bool surface_estimate(quadradic_surface_t &surface, float threshold);
    void planning_on_manifold_init();
    void surface_coordinate_to_traj_on_manifold(quadradic_surface_t &surface, traj_t &traj);
    void surface_differential(quadradic_surface_t &surface, Vector2f &coordinate_pos);
    float surface_coordinate_to_hight(float x, float y, Matrix<float, 6, 1> &para);
    void surface_search_box_set(quadradic_surface_t &surface);
    void set_motion_manifold_msg();
public:
    bool                                    traj_request;
    traj_t                                  traj;
    plane_t                                 tangent_plane;
    quadradic_surface_t                     surface;
    ekf_data_t                              ekf_data;
    std_msgs::Float32MultiArray             motion_manifold_msg;
    nav_msgs::Path                          motion_path;
    geometry_msgs::PoseStamped              motion_pose;

    void read_traj_in_EuclideanSpace(traj_t &traj, const traj_data_t &traj_data);
    void planning_on_manifold_main(KD_FOREST &ikdforest);
};

void traj_request_cb(const std_msgs::Bool::ConstPtr &msg);
#endif