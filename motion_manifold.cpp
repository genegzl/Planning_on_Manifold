#include "motion_manifold.hpp"
#include "traj_data.hpp"


extern const traj_data_t traj_data;


void traj_request_cb(const std_msgs::Bool::ConstPtr &msg)
{
    traj_request = *msg;
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

    traj.coordinate.yaw_rate = traj_data[traj.read_cnt][11];
    traj.coordinate.yaw = traj_data[traj.read_cnt][12];

    traj.coordinate.pos(0) += traj.pos_init(0);
    traj.coordinate.pos(1) += traj.pos_init(1);

    switch (traj.mode)
    {
        case TRAJ_REPEAT:
            traj.read_cnt ++;
            if (traj.read_cnt > TRAJ_DATA_ROW)
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
convert to: a1/a6*x^2 + a2/a6*y^2 + a3/a6*xy + a4/a6*x + a5/a6*y - 1/a6*z = -1
solve: A0*x0 = b0
where A0_i = [x_i^2, y_i^2, x_i*y_i, x_i, y_i, z_i], x0 = [a1/a6, a2/a6, a3/a6, a4/a6, a5/a6, -1/a6]^T, b0 = [-1, ..., -1]^T
*/
bool MANIFOLD::surface_estimate(quadradic_surface_t &surface, float threshold, int point_num)
{
    uint8_t j = 0;
    bool surfaceValid = true;
    float a6, esti_err;
    Matrix<float, 6,1> para_candidate;

    surface.lsm_b.setOnes();
    surface.lsm_b *= -1.0f;

    for (j = 0; j < point_num; j++)
    {
        surface.lsm_A(j,0) = surface.neighbor_ponits[j].x * surface.neighbor_ponits[j].x;
        surface.lsm_A(j,1) = surface.neighbor_ponits[j].y * surface.neighbor_ponits[j].y;
        surface.lsm_A(j,2) = surface.neighbor_ponits[j].x * surface.neighbor_ponits[j].y;
        surface.lsm_A(j,3) = surface.neighbor_ponits[j].x;
        surface.lsm_A(j,4) = surface.neighbor_ponits[j].y;
        surface.lsm_A(j,5) = surface.neighbor_ponits[j].z;
    }
    para_candidate = surface.lsm_A.colPivHouseholderQr().solve(surface.lsm_b);
    
    esti_err = (surface.lsm_A * para_candidate - surface.lsm_b).norm();

    if (esti_err > threshold)
    {
        surfaceValid = false;
    }
    else
    { 
        surface.para(5) = -1.0f / para_candidate(5);
        surface.para.block<5,1>(0,0) = para_candidate.block<5,1>(0,0) * surface.para(5);
    }

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

    /* 3. accleration */
    traj.acc.block<2,1>(0,0) = cdn->acc;
    traj.acc(2) = surface.p_F_x * cdn->acc(0) + surface.p_F_y * cdn->acc(1);

    /* 4. attitude */
    surface.normal_vector << surface.p_F_x, surface.p_F_y, -1.0f;
    surface.normal_vector.normalize();
    Bx << 1.0f, 0.0f, surface.p_F_x;
    Bx.normalize();
    By << 0.0f, 1.0f, surface.p_F_y;
    By.normalize();
    traj.R.col(0) = Bx;
    traj.R.col(1) = By;
    traj.R.col(2) = surface.normal_vector;
    traj.q = traj.R;
    traj.yaw = cdn->yaw;
    traj.omega(2) = cdn->yaw_rate;

    cout << traj.R << endl;
}

/* comment
F'_x: 2*a1*x + a3*y + a4
F''_xx: 2*a1
F'_y: 2*a2*y + a3*x + a5
F''_yy: 2*a2
normal vector:  [-F'_x, -F'_y, 1] / norm(.)
*/
void MANIFOLD::surface_differential(quadradic_surface_t &surface, Vector2f &coordinate_pos)
{
    surface.p_F_x = 2.0f*surface.para(0)*coordinate_pos(0) + surface.para(2)*coordinate_pos(1) + surface.para(3);
    surface.p_F_xx = 2.0f*surface.para(0);
    surface.p_F_y = 2.0f*surface.para(1)*coordinate_pos(1) + surface.para(2)*coordinate_pos(0) + surface.para(4);
    surface.p_F_yy = 2.0f*surface.para(1);
}

void MANIFOLD::planning_on_manifold_init()
{
    this->traj.init_en = true;
    this->traj.mode = TRAJ_REPEAT;
    this->traj.read_cnt = 0;
    this->traj.pos_init = this->ekf_data.pos;

    this->surface.search_center.x = this->traj.pos_init(0);
    this->surface.search_center.y = this->traj.pos_init(1);
    this->surface.search_center.z = this->traj.pos_init(2);
    this->init_en = true;
}

float MANIFOLD::surface_coordinate_to_hight(float x, float y, Matrix<float, 6, 1> &para)
{
    float z;
    z = para(0)*x*x + para(1)*y*y + para(2)*x*y + para(3)*x + para(4)*y + para(5);
    return z;
}

int db_cnt = 0;
void MANIFOLD::planning_on_manifold_main(KD_TREE &ikdtree)
{
    if (this->init_en != true)
    {
        this->planning_on_manifold_init();
    }

    this->read_traj_in_EuclideanSpace(this->traj, traj_data);
    ikdtree.Nearest_Search(this->surface.search_center, ESTI_SURFACE_PONIT_NUM, this->surface.neighbor_ponits, this->surface.Distance2Centrois);
    
    this->surfaceValid = this->surface_estimate(this->surface, SURFACE_FITTING_THRESHOLE, ESTI_SURFACE_PONIT_NUM);
    
    if (this->surfaceValid)
    {
        this->surface_differential(this->surface, this->traj.coordinate.pos);
        this->surface_coordinate_to_traj_on_manifold(this->surface, this->traj);
    }
    
    cout << "\n" << "surface para: " << this->surface.para.transpose() << endl;

    /* update next surface search center */
    this->surface.search_center.x = traj_data[this->traj.read_cnt][0];
    this->surface.search_center.y = traj_data[this->traj.read_cnt][1];
    this->surface.search_center.z = this->surface_coordinate_to_hight(this->surface.search_center.x, this->surface.search_center.y, this->surface.para);
    
    traj_data[this->traj.read_cnt][0];


    // std::cout << "TRAJ POSITION: " << this->traj.pos.transpose() << endl;
    // std::cout << "EKF POSITION: " << this->ekf_data.pos.transpose() << endl;
}