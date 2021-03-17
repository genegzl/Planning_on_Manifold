#ifndef  MOTION_MANIFOLD_HPP
#define  MOTION_MANIFOLD_HPP

#include <Eigen/Dense>
#include <opencv2/core.hpp>
#include <opencv2/core/eigen.hpp>
#include <ikd-Tree/ikd_Tree.h>

// using namespace std;
using namespace Eigen;


#define TRAJ_DATA_ROW 2001
#define TRAJ_DATA_COL 18
#define ESTI_PLANE_PONIT_NUM 5
#define ESTI_SURFACE_PONIT_NUM 125


typedef float traj_data_t[TRAJ_DATA_ROW][TRAJ_DATA_COL];

struct traj_t
{
    bool                                               init_en;
    // traj_mode_t                                        mode;
    // uint8_t                                            platform;
    // traj_state_t                                       state;
    uint32_t                                           wait_cnt;
    uint32_t                                           read_cnt;

    Vector3f                                           pos;                     // local position,      unit: m
    Vector3f                                           vel;                     // local velocity,      unit: m/s
    Vector3f                                           acc;                     // local acceleration,  unit: m/s^2
    Vector3f                                           omega;                   // body angular rate,   unit: rad/s
    float                                              yaw;                     // Euler yaw,           unit: rad
    float                                              a_T;                     // thrust acceleration, unit: N/kg or m/s^2
    Quaternionf                                        q;                       // quaternion,          order: wxyz

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
     Matrix<float, ESTI_PLANE_PONIT_NUM, 3>     lsm_A;
     Matrix<float, ESTI_PLANE_PONIT_NUM, 1>     lsm_b;
 };

 struct quadradic_surface_t
 {
     Matrix<float, 6, 1>                            para;
     PointVector                                    neighbor_ponits;
     Vector3f                                       centroid;
     Vector3f                                       normal_vector;
     Vector3f                                       x_proj;
     Vector3f                                       y_proj;
     Matrix<float, ESTI_SURFACE_PONIT_NUM, 6>       lsm_A;
     Matrix<float, ESTI_SURFACE_PONIT_NUM, 1>       lsm_b;
 };

 struct traj_manifold_t
 {
    traj_t                          traj;
    plane_t                         tangent_plane;
 };

class Manifold
{
private:
    bool                            init_en;
    void read_traj_in_EuclideanSpace(traj_t &traj, traj_data_t &traj_data);
    bool plane_estimate(plane_t &plane, float threshold, int point_num);
    void get_centroid(Vector3f &centroid, PointVector &point, int point_num);
    void point_prejoct_on_plane(Vector3f &p_out, Vector3f &p_in, plane_t &plane);
    bool surface_estimate(quadradic_surface_t &surface, float threshold, int point_num);
public:
    traj_t                          traj;
    plane_t                         tangent_plane;
    quadradic_surface_t             surface;
    void planning_on_manifold_main(Manifold &manifold);
};


#endif