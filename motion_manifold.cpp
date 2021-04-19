#include "motion_manifold.hpp"
#include "traj_data.hpp"

extern const traj_data_t traj_data;

float gene_math_sys_run_1st(float input, sys_1st_t &sys)
{
    if (sys.cfg_en == true)
    {
        float output = sys.num[0] * input + sys.num[1] * sys.x - sys.den[1] * sys.y;

        sys.x = input;
        sys.y = output;

        return output;
    }
    else
    {
        return input;
    } 
}

void math_butterworth_LPF_1st_config(uint8_t order, float Fc, float Ts, sys_1st_t &sys)
{
    float wc, wa, a;
    wc = PI_2 * Fc;
    wa = 2.0f / Ts * tanf(wc*Ts/2.0f);
    
    sys.num[0] = wa * Ts / (2.0f + wa * Ts);
    sys.num[1] = sys.num[0];
    sys.den[0] = 1.0f;
    sys.den[1] = (wa * Ts - 2.0f) / (2.0f + wa * Ts);
    sys.cfg_en = true;
}

/*
data: pos(0:2), vb_x, ab_x, yaw, yaw_rate, p_F_x, p_F_y, p_F_xx, p_F_yy
*/
void MANIFOLD::set_motion_manifold_msg()
{
    this->car_traj_msg.position.x = this->traj.pos(0);
    this->car_traj_msg.position.y = this->traj.pos(1);
    this->car_traj_msg.position.z = this->traj.pos(2);
    this->car_traj_msg.velocity_body.x = this->traj.vb;
    this->car_traj_msg.acceleration_body.x = this->traj.ab;
    this->car_traj_msg.yaw = this->traj.yaw;
    this->car_traj_msg.yaw_rate = this->traj.omega(2);
    this->car_traj_msg.p_F_x = this->surface.p_F_x;
    this->car_traj_msg.p_F_xx = this->surface.p_F_xx;
    this->car_traj_msg.p_F_y = this->surface.p_F_y;
    this->car_traj_msg.p_F_yy = this->surface.p_F_yy;
    this->car_traj_msg.p_F_xy = this->surface.p_F_xy;
    this->car_traj_msg.SurfaceValid = this->surface.surfaceValid;
}

uint32_t find_min_index(vector<float> &squence)
{
    uint32_t length, j, min_index;
    length = squence.size();
    min_index = 0;

    for (j = 1; j < length; j++)
    {
        if (squence[j] < squence[min_index])
        {
            min_index = j;
        }
    }

    return min_index;
}


void MANIFOLD::read_traj_in_EuclideanSpace(traj_t &traj, const traj_data_t &traj_data)
{
    
    if (traj.init_en != true)
    {
        return;
    }

    switch (traj.mode)
    {
        case TRAJ_REPEAT:
        {
            traj.read_cnt += traj.read_step;
            if (traj.read_cnt >= TRAJ_DATA_ROW)
            {
                traj.read_cnt = TRAJ_REPEAT_TICK;
            }
        }
        break;

        case TRAJ_STAY:
        {
            if (traj.read_cnt < TRAJ_DATA_ROW-1-traj.read_step)
            {
                traj.read_cnt += traj.read_step;
            }
            else
            {
                traj.read_cnt = TRAJ_DATA_ROW - 1;
            }
        }
        break;

        case TRAJ_MATCH:
        {
            uint32_t j, search_index;
            float px, diff;
            vector<float> px_sequence;

            if (traj.read_cnt < (TRAJ_CONSTANT_SPEED_START_TICK - traj.read_step) 
                || traj.read_cnt > (TRAJ_CONSTANT_SPEED_END_TICK - traj.read_step))
            {
                traj.constant_speed_enable = false;
                traj.read_cnt += traj.read_step;
            }
            else
            {
                traj.constant_speed_enable = true;

                px = traj.coordinate.pos(0) + traj.vel(0) * TRAJ_TS + 0.5f * traj.acc(0) * TRAJ_TS * TRAJ_TS - traj.pos_init(0);               
                for (j = 1; j <= TRAJ_MATCH_RANGE; j++)
                {
                    diff = fabs(traj_data[traj.read_cnt + j][0] - px);
                    px_sequence.push_back(diff);
                }
                search_index = find_min_index(px_sequence);

                // cout << "px: " << px << endl;
                // cout << "read_cnt: " << traj.read_cnt << endl;
                // cout << "vector size: " << px_sequence.size() << endl;
                // cout << "search index: " << search_index << endl;
                // cout << "search_list: " << endl;
                // for (j = 1; j <= px_sequence.size(); j++)
                // {
                //     cout << j << ", " << traj_data[traj.read_cnt + j][0] << ", " << px_sequence[j-1] << endl;
                // }

                traj.read_cnt += (search_index + 1);
                px_sequence.clear();
            }

            if (traj.read_cnt > TRAJ_DATA_ROW-1)
            {
                traj.read_cnt = TRAJ_DATA_ROW - 1;
            }

            
        }   
        break;

        default:
        break;
    }

    memcpy(traj.coordinate.pos.data(),   &traj_data[traj.read_cnt][0], 2*sizeof(float));
    memcpy(traj.coordinate.vel.data(),   &traj_data[traj.read_cnt][2], 2*sizeof(float));
    memcpy(traj.coordinate.acc.data(),   &traj_data[traj.read_cnt][4], 2*sizeof(float));

    traj.coordinate.vb = traj_data[traj.read_cnt][6];
    traj.coordinate.ab = traj_data[traj.read_cnt][7];
    traj.coordinate.yaw_rate = traj_data[traj.read_cnt][8];
    traj.coordinate.yaw = traj_data[traj.read_cnt][9];
    traj.coordinate.kp = traj_data[traj.read_cnt][10];

    traj.coordinate.pos(0) += traj.pos_init(0);
    traj.coordinate.pos(1) += traj.pos_init(1);

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
    bool surfaceValid = false;
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
        surface.lsm_A(j,0) = surface.neighbor_ponits[j].point.x * surface.neighbor_ponits[j].point.x;
        surface.lsm_A(j,1) = surface.neighbor_ponits[j].point.y * surface.neighbor_ponits[j].point.y;
        surface.lsm_A(j,2) = surface.neighbor_ponits[j].point.x * surface.neighbor_ponits[j].point.y;
        surface.lsm_A(j,3) = surface.neighbor_ponits[j].point.x;
        surface.lsm_A(j,4) = surface.neighbor_ponits[j].point.y;
        surface.lsm_A(j,5) = 1.0f;
        surface.lsm_b(j) = surface.neighbor_ponits[j].point.z;
    }

    para_candidate = surface.lsm_A.colPivHouseholderQr().solve(surface.lsm_b);

    surface.esti_err = (surface.lsm_A * para_candidate - surface.lsm_b).norm() / point_num;

    if (surface.esti_err > threshold)
    {
        surfaceValid = false;
    }
    else
    { 
        for (int j = 0; j < 6; j++)
        {
            surface.para[j] = gene_math_sys_run_1st(para_candidate[j], surface.fltr[j]);
        }
        surface.para = para_candidate;
    }

    // cout << "esti error: " << surface.esti_err << endl;
    // cout << "esti result: " << surfaceValid << endl;

    surfaceValid = true; // false; //
    return surfaceValid;
}

void MANIFOLD::surface_coordinate_to_traj_on_manifold(quadradic_surface_t &surface, traj_t &traj)
{
    coordinate_t *cdn = &traj.coordinate;
    Matrix<float, 6,1> para = surface.para;
    float v_norm, w, tmp1, tmp2, tmp0;

    if (surface.surfaceValid != true)
    {
        traj.pos.block<2,1>(0,0) = cdn->pos;
        traj.vel.block<2,1>(0,0) = cdn->vel;
        traj.vb = cdn->vb;
        traj.acc.block<2,1>(0,0) = cdn->acc;
        traj.ab = cdn->ab;
        traj.yaw = cdn->yaw;
        traj.omega(2) = cdn->yaw_rate;
        return;
    }
    else
    {
        if (traj.constant_speed_enable != true)
        {
            /* 1. postion */
            traj.pos.block<2,1>(0,0) = cdn->pos;
            traj.pos(2) = surface_coordinate_to_hight(traj.pos(0), traj.pos(1), surface.para);

            /* 2. velocity */
            traj.vel.block<2,1>(0,0) = cdn->vel;
            traj.vel(2) = surface.p_F_x * cdn->vel(0) + surface.p_F_y * cdn->vel(1);
            traj.vb = traj.vel.norm();
            if (cdn->vb < 0)
            {
                traj.vb *= -1.0f;
            }

            /* 3. accleration */
            traj.acc.block<2,1>(0,0) = cdn->acc;
            traj.acc(2) = surface.p_F_x * cdn->acc(0) + surface.p_F_y * cdn->acc(1);
            traj.ab = traj.acc.norm();
            if (cdn->ab < 0)
            {
                traj.ab *= -1.0f;
            }

            /* 4. attitude */
            traj.yaw = cdn->yaw;
            traj.omega(2) = cdn->yaw_rate;
        }
        else
        {
            /* 1. postion */
            traj.pos.block<2,1>(0,0) = cdn->pos;
            traj.pos(2) = surface_coordinate_to_hight(traj.pos(0), traj.pos(1), surface.para);

            /* 2. velocity */
            traj.vel.block<2,1>(0,0) = cdn->vel;
            traj.vel(2) = surface.p_F_x * cdn->vel(0) + surface.p_F_y * cdn->vel(1);
            v_norm = traj.vel.norm();
            traj.vel *= (CONSTANT_SPEED / v_norm);
            traj.vb = CONSTANT_SPEED;
            if (cdn->vb < 0)
            {
                traj.vb *= -1.0f;
            }

            /* 3. accleration */  
            /* constrains: (1). acc is vertical to vel (2). acc satisfies the surface equation (3). vel and acc on xy-plane satisfies the curvature on xy-plane */
            w = surface.p_F_xx * traj.vel(0) * traj.vel(0) + 2.0f * surface.p_F_xy * traj.vel(0) * traj.vel(1) + surface.p_F_yy * traj.vel(1) * traj.vel(1);
            tmp0 = traj.vel(0) * traj.vel(0) + traj.vel(1) * traj.vel(1);
            if (fabs(traj.vel(0)) > EPSILON)
            {
                tmp1 = cdn->kp * powf(tmp0, 1.5f);
                tmp2 = traj.vel(1) + traj.vel(2) * surface.p_F_y;
                traj.acc(0) = (-tmp2*tmp1/traj.vel(0) - traj.vel(2)*w) / (traj.vel(0) + traj.vel(2)*surface.p_F_x + tmp2*traj.vel(1)/traj.vel(0));
                traj.acc(1) = (tmp1 + traj.vel(1)*traj.acc(0)) / traj.vel(0);
            }
            else
            {
                traj.acc(0) = - cdn->kp * traj.vel(1) * traj.vel(1);
                traj.acc(1) = - traj.vel(2) * (w + surface.p_F_x*traj.acc(0)) / (traj.vel(1) * surface.p_F_y);
            }
            traj.acc(2) = w + surface.p_F_x * traj.acc(0) + surface.p_F_y * traj.acc(1);
            traj.ab = traj.acc.norm();
            // todo: when ab < 0 ?

            // cout << "surface partical: " << surface.p_F_x << ", " << surface.p_F_y << ", " << surface.p_F_xx << ", " << surface.p_F_xy << ", " << surface.p_F_yy << endl;
            // cout << "curvature on plane: " << cdn->kp << endl;
            // cout << "coordinate vel: " << cdn->vel.transpose() << endl;
            // cout << "traj vel: " << traj.vel.transpose() << endl;
            // cout << "traj acc: " << traj.acc.transpose() << endl;
            // cout << "w: " << w << endl;
            // cout << "tmp0: " << tmp0 << endl;
            // cout << "tmp1: " << tmp1 << endl;
            // cout << "tmp2: " << tmp2 << endl; 

            /* 4. attitude */
            traj.yaw = atan2(traj.vel(1), traj.vel(0)); // CONSTANT_SPEED != 0

            /* 5. yaw rate */
            traj.omega(2) = (traj.vel(0) * traj.acc(1) - traj.acc(0) * traj.vel(1)) / tmp0;
        }
        
    } 
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
    if (surface.surfaceValid == true)
    {
        surface.search_box.vertex_min[1] = surface.search_center(1) - SURFACE_SEARCH_RANGE_RIGHT;
        surface.search_box.vertex_min[2] = surface.search_center(2) - SURFACE_SEARCH_RANGE_DOWN;
    }
    else
    {
        surface.search_box.vertex_min[1] = surface.search_center(1) - 2.0f * SURFACE_SEARCH_RANGE_RIGHT;
        surface.search_box.vertex_min[2] = surface.search_center(2) - 2.0f * SURFACE_SEARCH_RANGE_DOWN;
    }
    
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
    this->traj.mode = TRAJ_MATCH;
    this->traj.read_step = TRAJ_DATA_FS / TRAJ_OUPUT_FS;
    this->traj.read_cnt = 0;
    this->traj.pos_init = this->ekf_data.pos;
    this->traj.pos_init.setZero();

    this->surface.search_center = this->traj.pos_init;
    this->surface.search_center(2) -= LIDAR_MOUNT_HIGHT;

    for (int j = 0; j < 6; j++)
    {
        math_butterworth_LPF_1st_config(1, 0.5f, 0.02f, this->surface.fltr[j]);
    }

    this->init_en = true;
}

int db_cnt = 0;
float p_integ_db, v_integ_db;
int test_axis = 0;
void MANIFOLD::planning_on_manifold_main(KD_FOREST &ikdforest)
{
    if (this->init_en != true)
    {
        this->planning_on_manifold_init();
    }

    this->read_traj_in_EuclideanSpace(this->traj, traj_data);

    this->surface_search_box_set(this->surface);

    this->surface.neighbor_ponits.clear();
    ikdforest.Box_Search(this->surface.search_box, this->surface.neighbor_ponits);

    this->surface.surfaceValid = this->surface_estimate(this->surface, SURFACE_FITTING_THRESHOLE);

    if (this->surface.surfaceValid)
    {
        this->surface_differential(this->surface, this->traj.coordinate.pos);
    }

    this->surface_coordinate_to_traj_on_manifold(this->surface, this->traj);
    this->set_motion_manifold_msg();

    if ( db_cnt != 0) //
    {
      
        // cout << "neighbor points: " << endl;
        // int point_num;
        // if (surface.neighbor_point_num > SURFACE_ESTI_PONIT_NUM_MAX)
        // {
        //     point_num = SURFACE_ESTI_PONIT_NUM_MAX;
        // }
        // else
        // {
        //     point_num = surface.neighbor_point_num;
        // }
        
        // for (int j = 0; j < point_num; j++)
        // {
        //     cout << surface.neighbor_ponits[j].x << ", " << surface.neighbor_ponits[j].y << ", " << surface.neighbor_ponits[j].z << endl;
        // }
        
        cout << "surface valid: " << surface.surfaceValid << endl;
        cout << "ekf position: " << ekf_data.pos.transpose() << endl;
        cout << "search center: " << surface.search_center.transpose() << endl;
        // cout << "search box min: " << surface.search_box.vertex_min[0] << ", " << surface.search_box.vertex_min[1] << ", " << surface.search_box.vertex_min[2] << endl;
        // cout << "search box max: " << surface.search_box.vertex_max[0] << ", " << surface.search_box.vertex_max[1] << ", " << surface.search_box.vertex_max[2] << endl;
        cout << "neighbor point num: " << surface.neighbor_point_num << endl;
        // cout << "traj coordinate: " << traj.coordinate.pos.transpose() << endl;
        cout << "esti error: " << surface.esti_err << endl;
        cout << "fltr num" << surface.fltr->num[0] << ", " << surface.fltr->num[1] << endl;
        cout << "fltr den" << surface.fltr->den[0] << ", " << surface.fltr->den[1] << endl;
        // cout << traj.pos(0) << ", " << traj.pos(1) << ", " << traj.pos(2) << ", " << traj.vb << ", " << traj.yaw << ", " << traj.omega(2) << ", " << surface.p_F_x << ", " << surface.p_F_y << ", " << surface.p_F_xx << ", " << surface.p_F_xy << ", " << surface.p_F_yy << endl;

        // cout << "traj pos: " << traj.pos.transpose() << endl;
        // cout << "traj init pos: " << traj.pos_init.transpose() << endl;
        // cout << "traj data:pos: " << traj_data[this->traj.read_cnt-1][0] << ", " << traj_data[this->traj.read_cnt-1][1] << endl; 

    }

    db_cnt ++;
    
    if (db_cnt == 1)
    {
        v_integ_db = traj.vel(test_axis);
        p_integ_db = traj.pos(test_axis);
    }
    else
    {
        p_integ_db += traj.vel(test_axis) * TRAJ_TS;
        // p_integ_db += v_integ_db * TRAJ_TS + 0.0f * traj.acc(test_axis) * TRAJ_TS * TRAJ_TS;
        v_integ_db += traj.acc(test_axis) * TRAJ_TS;
    }

    float v_norm_db = traj.vel.norm();
    
    db_p.push_back(traj.pos(test_axis));  
    db_p_integ.push_back(p_integ_db);
    db_v.push_back(traj.vel(test_axis));  
    db_v_integ.push_back(v_integ_db);
    db_v_norm.push_back(v_norm_db);
    db_time_stamp.push_back(db_cnt * TRAJ_TS);

    /* update next surface search center */
    this->surface.search_center(0) = traj_data[this->traj.read_cnt][0] + traj.pos_init(0);
    this->surface.search_center(1) = traj_data[this->traj.read_cnt][1] + traj.pos_init(1);
    if (surface.surfaceValid == true)
    {
        this->surface.search_center(2) = this->surface_coordinate_to_hight(this->surface.search_center(0), this->surface.search_center(1), this->surface.para);
    }
    

    // std::cout << "TRAJ POSITION: " << this->traj.pos.transpose() << endl;
    // std::cout << "EKF POSITION: " << this->ekf_data.pos.transpose() << endl;
}