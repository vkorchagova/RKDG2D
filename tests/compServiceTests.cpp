/// Test for compService

#include <vector>
#include <memory>
#include "compService.h"
#include "numvector.h"
#include "Point.h"
#include "Params.h"
#include "gtest/gtest.h"

namespace {

::testing::AssertionResult equalNumvectors(const numvector<double, dimPh> v1, const numvector<double, dimPh> v2)
{
    for (size_t i = 0; i < dimPh; ++i)
    {
        if ( v1[i] != v2[i]){
            return ::testing::AssertionFailure() << "array[" << i
                << "] (" << v1[i] << ") != expected[" << i
                << "] (" << v2[i] << ")";
        }
    }

    return ::testing::AssertionSuccess();
}

TEST(CompService, rotation) {

    numvector<double, dimPh> sol = {1.0, 1.0, 0.0, 0.0, 1.0};
    Point n ({0.0, 1.0});

    numvector<double, dimPh> result = rotate(sol, n);
    numvector<double, dimPh> ref = {1.0, 0.0, -1.0, 0.0, 1.0};

    EXPECT_TRUE(equalNumvectors(result, ref));
}

TEST(CompService, rotationInverse) {

    numvector<double, dimPh> sol = {1.0, 1.0, 0.0, 0.0, 1.0};
    Point n ({0.0, 1.0});

    numvector<double, dimPh> result = inverseRotate(sol, n);
    numvector<double, dimPh> ref = {1.0, 0.0, 1.0, 0.0, 1.0};

    EXPECT_TRUE(equalNumvectors(result, ref));
}

TEST(CompService, integrateScalarFunction) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.0}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::function<double(const Point &)> f = [&](const Point& p){ return 1.0; };

    double res = integrate(e, f);

    EXPECT_FLOAT_EQ(res, e.getLength());
}

TEST(CompService, integrateScalarFunctionWithGivenValues) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.0}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::vector<double> f (e.nGP);
    for (int i = 0; i < e.nGP; ++i)
        f[i] = 1.0;

    double res = integrate(e, f);

    EXPECT_FLOAT_EQ(res, e.getLength());
}

}  // namespace