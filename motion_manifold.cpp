#include "motion_manifold.hpp"
#include "traj_data.hpp"


extern const traj_data_t traj_data;

void Manifold::read_traj_in_EuclideanSpace(traj_t &traj, traj_data_t &traj_data)
{
    memcpy(traj.pos.data(),   &traj_data[traj.read_cnt][0], 3*sizeof(float));
    memcpy(traj.vel.data(),   &traj_data[traj.read_cnt][3], 3*sizeof(float));
    memcpy(traj.acc.data(),   &traj_data[traj.read_cnt][6], 3*sizeof(float));
    memcpy(traj.omega.data(), &traj_data[traj.read_cnt][9], 3*sizeof(float));

    traj.yaw = traj_data[traj.read_cnt][12];
    // traj.a_T = traj_data[traj.read_cnt][13];

    traj.q.w() = traj_data[traj.read_cnt][14];
    traj.q.x() = traj_data[traj.read_cnt][15];
    traj.q.y() = traj_data[traj.read_cnt][16];
    traj.q.z() = traj_data[traj.read_cnt][17];
}

/* comment
plane equation: Ax + By + Cz + D = 0
convert to: A/D*x + B/D*y + C/D*z = -1
solve: A0*x0 = b0
where A0_i = [x_i, y_i, z_i], x0 = [A/D, B/D, C/D]^T, b0 = [-1, ..., -1]^T
normvec:  normalized x0
*/
bool Manifold::plane_estimate(plane_t &plane, float threshold, int point_num)
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

void Manifold::get_centroid(Vector3f &centroid, PointVector &point, int point_num)
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

void Manifold::point_prejoct_on_plane(Vector3f &p_out, Vector3f &p_in, plane_t &plane)
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
F'_x: 2*a1*x + a3*y + a4
F''_xx: 2*a1
F'_y: 2*a2*y + a3*x + a5
F''_yy: 2*a2
normal vector:  [-F'_x, -F'_y, 1] / norm(.)
*/
bool Manifold::surface_estimate(quadradic_surface_t &surface, float threshold, int point_num)
{
    uint8_t j = 0;
    bool surfaceValid = true;
    float a6, esti_err;

    surface.lsm_b.setOnes();
    surface.lsm_b *= -1.0f;

    for (j = 0; j < point_num; j++)
    {
        surface.lsm_A(j,0) = surface.neighbor_ponits[j].x;
        surface.lsm_A(j,1) = surface.neighbor_ponits[j].y;
        surface.lsm_A(j,2) = surface.neighbor_ponits[j].z;
    }
    surface.para = surface.lsm_A.colPivHouseholderQr().solve(surface.lsm_b);
    
    esti_err = (surface.lsm_A * surface.para - surface.lsm_b).norm();

    if (esti_err > threshold)
    {
        surfaceValid = false;
    }

    return surfaceValid;
}


void Manifold::planning_on_manifold_main(Manifold &manifold)
{
    if (manifold.init_en != true)
    {
        manifold.init_en = true;
    }
}