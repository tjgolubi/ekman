#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/srs/epsg.hpp>
#include <boost/geometry/srs/esri.hpp>
#include <boost/geometry/srs/iau2000.hpp>
#include <boost/geometry/srs/projection.hpp>

#include <iostream>
#include <iomanip>

int main(int, char*[]) {
    using namespace boost::geometry;
    using namespace boost::geometry::model;
    using namespace boost::geometry::srs;

    using point_ll =  point<double, 2, cs::geographic<degree>>;
    using point_xy =  point<double, 2, cs::cartesian>;

    using namespace boost::geometry::srs::dpar;

    const point_ll pt_ll0(-91.661078, 41.907199);
    const point_ll pt_ll (-91.655476, 41.903958);
    const point_xy pt_xy0(0, 0);

    std::cout << " dist_ll=" << distance(pt_ll, pt_ll0) << std::endl;

    point_xy pt_xy(0, 0);
    point_ll pt_ll2(0, 0);

    projection<> prj = parameters<>(proj_aeqd)(ellps_wgs84)(lat_0, 41.907199)(lon_0, -91.661078)(x_0,0)(y_0,0)(units_m);

    prj.forward(pt_ll, pt_xy);

    std::cout << get<0>(pt_xy) << ' ' << get<1>(pt_xy) << " dist=" << distance(pt_xy, pt_xy0) << std::endl;

    prj.inverse(pt_xy, pt_ll2);

    std::cout << get<0>(pt_ll2) << ' ' << get<1>(pt_ll2) << " dist=" << distance(pt_ll2, pt_ll) << std::endl;

    return 0;
}
