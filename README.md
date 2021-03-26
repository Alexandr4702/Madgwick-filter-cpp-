U can use any vector for getting orientation such as acclereation, magnetick field, sun sensor vector.
Using example:
~~~
#include <iostream>
#include <fstream>
#include "Madgwick.h"

using namespace std;

const double degre_to_rad = M_PI / 180.0;
const double rad_to_degre = 180.0 / M_PI;

std::ostream & operator << (std::ostream & s, Quat q)
{
    s << q.w() << " " << q.x() << " " << q.y() << " " << q.z();
    return s;
}

int main()
{
    //uint16_t year, uint8_t month, uint8_t day, uint8_t hour, uint8_t minute, double second
    double start_time = 0;
    double current_time;
    double finish_time = 1000;
    double time_in_seconds = finish_time - start_time;
    double dt = 0.1;

    current_time = start_time;

    Madgwick_filter madgwick;
    Vec3 ref1(1,0,0);
    Vec3 ref2(0.5,0,0.5);
    ref2.normalize();
    Quat ref1_quat(0,ref1.x(), ref1.y(), ref1.z());
    Quat ref2_quat(0,ref2.x(), ref2.y(), ref2.z());
    Quat start_orienataion (1,0,0,0);
    Quat finish_orientation (0.924, 0, 0.383, 0);
    Quat current_orientaion = start_orienataion;

    Vec3 TrueOmega = acos(finish_orientation.w()) * 2 * finish_orientation.vec().normalized() / time_in_seconds;
    Vec3 meas1 = (current_orientaion * ref1_quat * current_orientaion.inverse()).vec();
    Vec3 meas2 = (current_orientaion * ref2_quat * current_orientaion.inverse()).vec();
    Vec3 omega_bias = Vec3(0, 0.1, 0.4) * degre_to_rad;
    Vec3 omega_meas = TrueOmega + omega_bias;
    uint64_t n = 0;
    while(current_time < finish_time)
    {
        madgwick.update(
                    ref1, ref2,
                    meas1, meas2,
                    omega_meas, dt);
//        cerr << madgwick.get_orientaion() << " madgwick " << n <<endl;
//        cerr << current_orientaion << " current " << n <<endl;
        cout << 2 * (madgwick.get_orientaion() - current_orientaion).norm() << " difference " << n << " n ";
        cout << madgwick.getOmega_bias().transpose() * rad_to_degre << " bias " << madgwick.get_orientaion() << " madgwick " << current_orientaion << " current " << endl;

        current_time += dt;
        Quat omega_quat;
        double omega_2 = (TrueOmega.transpose() * TrueOmega);
        double second_order = omega_2 * dt * dt / 8.0;
        omega_quat.w() = 1.0 - second_order;
        omega_quat.vec() = dt * 0.5 * TrueOmega;
        current_orientaion = current_orientaion * omega_quat;
        meas1 = (current_orientaion * ref1_quat * current_orientaion.inverse()).vec();
        meas2 = (current_orientaion * ref2_quat * current_orientaion.inverse()).vec();
        n++;
    }
}
~~~
