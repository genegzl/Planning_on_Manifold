#include "motion_manifold.hpp"
#include "traj_data.hpp"


extern const traj_data_t traj_data;


void traj_request_cb(const std_msgs::Bool::ConstPtr &msg)
{
    traj_request = *msg;
}

/*
data: pos(0:2), vb_x, ab_x, yaw, yaw_rate, p_F_x, p_F_y, p_F_xx, p_F_yy
*/
void MANIFOLD::set_motion_manifold_msg()
{
    this->motion_manifold_msg.data.clear();

    this->motion_manifold_msg.data.push_back(this->traj.pos(0));
    this->motion_manifold_msg.data.push_back(this->traj.pos(1));
    this->motion_manifold_msg.data.push_back(this->traj.pos(2));
    this->motion_manifold_msg.data.push_back(this->traj.vb);
    this->motion_manifold_msg.data.push_back(this->traj.ab);
    this->motion_manifold_msg.data.push_back(this->traj.yaw);
    this->motion_manifold_msg.data.push_back(this->traj.omega(2));
    this->motion_manifold_msg.data.push_back(this->surface.p_F_x);
    this->motion_manifold_msg.data.push_back(this->surface.p_F_y);
    this->motion_manifold_msg.data.push_back(this->surface.p_F_xx);
    this->motion_manifold_msg.data.push_back(this->surface.p_F_yy);
}

void MANIFOLD::read_traj_in_EuclideanSpace(traj_t &traj, const traj_data_t &traj_data)
{
    if (traj.init_en != true)
    {
        return;
    }

    memcpy(traj.coordinate.pos.data(),   &traj_data[traj.read_cnt][0], 2*sizeof(float));
    memcpy(traj.coordinate.vel.data(),   &traj_data[traj.read_cnt][3], 2*sizeof(float));
    memcpy(traj.coordinate.acc.data(),   &traj_data[traj.read_cnt][6], 2*sizeof(float));

    traj.coordinate.vb = traj_data[traj.read_cnt][9];
    traj.coordinate.ab = traj_data[traj.read_cnt][10];
    traj.coordinate.yaw_rate = traj_data[traj.read_cnt][11];
    traj.coordinate.yaw = traj_data[traj.read_cnt][12];

    traj.coordinate.pos(0) += traj.pos_init(0);
    traj.coordinate.pos(1) += traj.pos_init(1);

    switch (traj.mode)
    {
        case TRAJ_REPEAT:
            traj.read_cnt ++;
            if (traj.read_cnt >= TRAJ_DATA_ROW)
            {
                traj.read_cnt = TRAJ_REPEAT_TICK;
            }
            break;

        case TRAJ_STAY:
            if (traj.read_cnt < TRAJ_DATA_ROW)
            {
                traj.read_cnt ++;
            }
            break;

        default:
            break;
    }
}

/* comment
plane equation: Ax + By + Cz + D = 0
convert to: A/D*x + B/D*y + C/D*z = -1
solve: A0*x0 = b0
where A0_i = [x_i, y_i, z_i], x0 = [A/D, B/D, C/D]^T, b0 = [-1, ..., -1]^T
normvec:  normalized x0
*/
bool MANIFOLD::plane_estimate(plane_t &plane, float threshold, int point_num)
{
    uint8_t j = 0;
    bool planeValid = true;

    plane.lsm_b.setOnes();
    plane.lsm_b *= -1.0f;

    for (j = 0; j < point_num; j++)
    {
        plane.lsm_A(j,0) = plane.neighbor_ponits[j].x;
        plane.lsm_A(j,1) = plane.neighbor_ponits[j].y;
        plane.lsm_A(j,2) = plane.neighbor_ponits[j].z;
    }
    plane.normal_vector = plane.lsm_A.colPivHouseholderQr().solve(plane.lsm_b);
    
    for (j = 0; j < point_num; j++)
    {
        if (fabs(plane.normal_vector(0) * plane.neighbor_ponits[j].x + plane.normal_vector(1) * plane.neighbor_ponits[j].y + plane.normal_vector(2) * plane.neighbor_ponits[j].z + 1.0f) > threshold)
        {
            planeValid = false;
            return planeValid;
            break;
        }
    }

    plane.normal_vector.normalize();
    return planeValid;
}

void MANIFOLD::get_centroid(Vector3f &centroid, PointVector &point, int point_num)
{
    uint8_t j = 0;

    centroid.setZero();
    for (j = 0; j < point_num; j++)
    {
        centroid(0) += point[j].x;
        centroid(1) += point[j].y;
        centroid(2) += point[j].z;
    }
    centroid /= point_num;
}

void MANIFOLD::point_prejoct_on_plane(Vector3f &p_out, Vector3f &p_in, plane_t &plane)
{
    /* solve: dot((p_out - p_in) , (p_out - c))
              (px_out - px-in)/A = (py_out - py-in)/B = (pz_out - pz-in)/C
    */

    float AA = plane.normal_vector(0)*plane.normal_vector(0); // (B^2 + C^2) / A

    p_out(0) = plane.normal_vector(0) * (plane.normal_vector(0) * plane.centroid(0) + plane.normal_vector(1)*(plane.centroid(1)-p_in(1))
                + plane.normal_vector(2)*(plane.centroid(2)-p_in(2))) + (1.0f-AA)*p_in(0);
    p_out(1) = plane.normal_vector(1) * (p_out(0) - p_in(0)) / plane.normal_vector(0) + p_in(1);
    p_out(2) = plane.normal_vector(2) * (p_out(0) - p_in(0)) / plane.normal_vector(0) + p_in(2);
}

// void project_motion_to_manifold(Manifold &mm)
// {
//     Vector3f p = mm.traj.pos;
// }


/* comment
surface equation: z = F(x,y) = a1*x^2 + a2*y^2 + a3*xy + a4*x + a5*y + a6
solve: A0*x0 = b0
where A0_i = [x_i^2, y_i^2, x_i*y_i, x_i, y_i, 1], x0 = [a1, a2, a3, a4, a5, a6]^T, b0_i = z_i
*/
bool MANIFOLD::surface_estimate(quadradic_surface_t &surface, float threshold)
{
    uint8_t j = 0;
    uint32_t point_num;
    bool surfaceValid = true;
    Matrix<float, 6,1> para_candidate;

    surface.neighbor_point_num = surface.neighbor_ponits.size();

    if (surface.neighbor_point_num < SURFACE_ESTI_PONIT_NUM_MIN)
    {
        // cout << "The number of ground points is not enough!!!!!!!!!" << endl;
        return false;
    }
    else if (surface.neighbor_point_num > SURFACE_ESTI_PONIT_NUM_MAX)
    {
        point_num = SURFACE_ESTI_PONIT_NUM_MAX;
    }
    else
    {
        point_num = surface.neighbor_point_num;
    }

    surface.lsm_A.resize(point_num, 6);
    surface.lsm_b.resize(point_num, 1);

    for (j = 0; j < point_num; j++)
    {
        surface.lsm_A(j,0) = surface.neighbor_ponits[j].x * surface.neighbor_ponits[j].x;
        surface.lsm_A(j,1) = surface.neighbor_ponits[j].y * surface.neighbor_ponits[j].y;
        surface.lsm_A(j,2) = surface.neighbor_ponits[j].x * surface.neighbor_ponits[j].y;
        surface.lsm_A(j,3) = surface.neighbor_ponits[j].x;
        surface.lsm_A(j,4) = surface.neighbor_ponits[j].y;
        surface.lsm_A(j,5) = 1.0f;
        surface.lsm_b(j) = surface.neighbor_ponits[j].z;
    }

    para_candidate = surface.lsm_A.colPivHouseholderQr().solve(surface.lsm_b);

    surface.esti_err = (surface.lsm_A * para_candidate - surface.lsm_b).norm();

    if (surface.esti_err > threshold)
    {
        surfaceValid = false;
    }
    else
    { 
        surface.para = para_candidate;
    }

    // cout << "esti error: " << surface.esti_err << endl;
    // cout << "esti result: " << surfaceValid << endl;

    return surfaceValid;
}

void MANIFOLD::surface_coordinate_to_traj_on_manifold(quadradic_surface_t &surface, traj_t &traj)
{
    coordinate_t *cdn = &traj.coordinate;
    Matrix<float, 6,1> para = surface.para;
    Vector3f Bx, By;

    /* 1. postion */
    traj.pos.block<2,1>(0,0) = cdn->pos;
    traj.pos(2) = surface_coordinate_to_hight(traj.pos(0), traj.pos(1), surface.para);

    /* 2. velocity */
    traj.vel.block<2,1>(0,0) = cdn->vel;
    traj.vel(2) = surface.p_F_x * cdn->vel(0) + surface.p_F_y * cdn->vel(1);
    traj.vb = traj.vel.norm();
    if (traj.coordinate.vb < 0)
    {
        traj.vb *= -1.0f;
    }

    /* 3. accleration */
    traj.acc.block<2,1>(0,0) = cdn->acc;
    traj.acc(2) = surface.p_F_x * cdn->acc(0) + surface.p_F_y * cdn->acc(1);
    traj.ab = traj.coordinate.ab / fabs(traj.coordinate.ab) * traj.acc.norm();
    if (traj.coordinate.ab < 0)
    {
        traj.ab *= -1.0f;
    }

    /* 4. attitude */
    surface.normal_vector << -surface.p_F_x, -surface.p_F_y, 1.0f;
    surface.normal_vector.normalize();
    Bx << 1.0f, 0.0f, surface.p_F_x;
    Bx.normalize();
    By << 0.0f, 1.0f, surface.p_F_y;
    By.normalize();

    traj.yaw = cdn->yaw;
    traj.omega(2) = cdn->yaw_rate;
}

/* comment
F'_x: 2*a1*x + a3*y + a4
F''_xx: 2*a1
F'_y: 2*a2*y + a3*x + a5
F''_yy: 2*a2
F''_xy: a3
normal vector:  [-F'_x, -F'_y, 1] / norm(.)
*/
void MANIFOLD::surface_differential(quadradic_surface_t &surface, Vector2f &coordinate_pos)
{
    surface.p_F_x = 2.0f * surface.para(0)*coordinate_pos(0) + surface.para(2)*coordinate_pos(1) + surface.para(3);
    surface.p_F_xx = 2.0f * surface.para(0);
    surface.p_F_y = 2.0f * surface.para(1)*coordinate_pos(1) + surface.para(2)*coordinate_pos(0) + surface.para(4);
    surface.p_F_yy = 2.0f * surface.para(1);
    surface.p_F_xy = surface.para(2);
}

void MANIFOLD::surface_search_box_set(quadradic_surface_t &surface)
{
    surface.search_box.vertex_max[0] = surface.search_center(0) + SURFACE_SEARCH_RANGE_FORWARD;
    surface.search_box.vertex_max[1] = surface.search_center(1) + SURFACE_SEARCH_RANGE_LEFT;
    surface.search_box.vertex_max[2] = surface.search_center(2) + SURFACE_SEARCH_RANGE_UP;
    surface.search_box.vertex_min[0] = surface.search_center(0) - SURFACE_SEARCH_RANGE_BACK;
    surface.search_box.vertex_min[1] = surface.search_center(1) - SURFACE_SEARCH_RANGE_RIGHT;
    surface.search_box.vertex_min[2] = surface.search_center(2) - SURFACE_SEARCH_RANGE_DOWN;
}

float MANIFOLD::surface_coordinate_to_hight(float x, float y, Matrix<float, 6, 1> &para)
{
    float z;
    z = para(0)*x*x + para(1)*y*y + para(2)*x*y + para(3)*x + para(4)*y + para(5);
    return z;
}

void MANIFOLD::planning_on_manifold_init()
{
    this->traj.init_en = true;
    this->traj.mode = TRAJ_STAY;
    this->traj.read_cnt = 0;
    this->traj.pos_init = this->ekf_data.pos;

    this->traj.pos_init.setZero();

    this->surface.search_center = this->traj.pos_init;
    this->surface.search_center(2) -= LIDAR_MOUNT_HIGHT;

    this->init_en = true;
}

int db_cnt = 0;
void MANIFOLD::planning_on_manifold_main(KD_FOREST &ikdforest)
{
    if (this->init_en != true)
    {
        this->planning_on_manifold_init();
    }

    this->read_traj_in_EuclideanSpace(this->traj, traj_data);

    this->surface_search_box_set(this->surface);

    ikdforest.Box_Search(this->surface.search_box, this->surface.neighbor_ponits);

    this->surfaceValid = this->surface_estimate(this->surface, SURFACE_FITTING_THRESHOLE);

    if (this->surfaceValid)
    {
        this->surface_differential(this->surface, this->traj.coordinate.pos);
    }

    this->surface_coordinate_to_traj_on_manifold(this->surface, this->traj);
    this->set_motion_manifold_msg();

    if ( db_cnt != 0) //
    {
      
        cout << "neighbor points: " << endl;
        int point_num;
        if (surface.neighbor_point_num > SURFACE_ESTI_PONIT_NUM_MAX)
        {
            point_num = SURFACE_ESTI_PONIT_NUM_MAX;
        }
        else
        {
            point_num = surface.neighbor_point_num;
        }
        
        for (int j = 0; j < point_num; j++)
        {
            cout << surface.neighbor_ponits[j].x << ", " << surface.neighbor_ponits[j].y << ", " << surface.neighbor_ponits[j].z << endl;
        }
        
        
        cout << "ekf position: " << ekf_data.pos.transpose() << endl;
        cout << "search center: " << surface.search_center.transpose() << endl;
        cout << "search box min: " << surface.search_box.vertex_min[0] << ", " << surface.search_box.vertex_min[1] << ", " << surface.search_box.vertex_min[2] << endl;
        cout << "search box max: " << surface.search_box.vertex_max[0] << ", " << surface.search_box.vertex_max[1] << ", " << surface.search_box.vertex_max[2] << endl;
        cout << "neighbor point num: " << surface.neighbor_point_num << endl;
        // cout << "\n" << "surface para: " << this->surface.para.transpose() << endl;
        // cout << "\n" << "traj R: " << this->traj.R << endl;
        // cout << "\n" << "ekf R: " << this->ekf_data.R << endl;
        cout << "traj coordinate: " << traj.coordinate.pos.transpose() << endl;
        cout << "esti error: " << surface.esti_err << endl;
        // cout << "traj pos: " << traj.pos.transpose() << endl;
        // cout << "traj init pos: " << traj.pos_init.transpose() << endl;
        // cout << "traj data:pos: " << traj_data[this->traj.read_cnt-1][0] << ", " << traj_data[this->traj.read_cnt-1][1] << endl; 

    }

    db_cnt ++;
    
    // if (db_cnt < 0)
    // {
    //     cout << db_cnt << endl;
    // }
   

    /* update next surface search center */
    this->surface.search_center(0) = traj_data[this->traj.read_cnt][0] + traj.pos_init(0);
    this->surface.search_center(1) = traj_data[this->traj.read_cnt][1] + traj.pos_init(1);
    this->surface.search_center(2) = this->surface_coordinate_to_hight(this->surface.search_center(0), this->surface.search_center(1), this->surface.para);
    this->surface.neighbor_ponits.clear();

    // std::cout << "TRAJ POSITION: " << this->traj.pos.transpose() << endl;
    // std::cout << "EKF POSITION: " << this->ekf_data.pos.transpose() << endl;
}