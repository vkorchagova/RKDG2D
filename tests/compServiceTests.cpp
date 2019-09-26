/// Tests for compService

#include <vector>
#include <memory>
#include "compService.h"
#include "numvector.h"
#include "Point.h"
#include "Params.h"
#include "gtest/gtest.h"

namespace {

// check if two numvectors are equal
::testing::AssertionResult equalNumvectors(const numvector<double, dimPh> v1, const numvector<double, dimPh> v2)
{
    for (size_t i = 0; i < dimPh; ++i)
    {
        if ( fabs(v1[i] - v2[i]) > 1e-12){
            return ::testing::AssertionFailure() << "array[" << i
                << "] (" << v1[i] << ") != expected[" << i
                << "] (" << v2[i] << ")";
        }
    }

    return ::testing::AssertionSuccess();
}

// Test rotation
TEST(CompService, rotation) {

    numvector<double, dimPh> sol = {1.0, 1.0, 0.0, 0.0, 1.0};
    Point n ({0.0, 1.0});

    numvector<double, dimPh> result = rotate(sol, n);
    numvector<double, dimPh> ref = {1.0, 0.0, -1.0, 0.0, 1.0};

    EXPECT_TRUE(equalNumvectors(result, ref));
}

// Test inverse rotation
TEST(CompService, rotationInverse) {

    numvector<double, dimPh> sol = {1.0, 1.0, 0.0, 0.0, 1.0};
    Point n ({0.0, 1.0});

    numvector<double, dimPh> result = inverseRotate(sol, n);
    numvector<double, dimPh> ref = {1.0, 0.0, 1.0, 0.0, 1.0};

    EXPECT_TRUE(equalNumvectors(result, ref));
}

// Test of integration of the scalar function through 30 deg rotated edge 
// Input data:
// - edge (0;0) - (0.5*sqrt(3), 0.5)
// - test function (returns 1.0 in all points) defined as std::function
// Success: result should be equal to the edge length
TEST(CompService, integrateScalarFunction1D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.5}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::function<double(const Point &)> f = [&](const Point& p){ return 1.0; };

    double res = integrate(e, f);

    EXPECT_FLOAT_EQ(res, e.getLength());
}

// Test of integration of the scalar function through 30 deg rotated edge 
// Input data:
// - edge (0;0) - (0.5*sqrt(3), 0.5)
// - test function (returns 1.0 in all points) values in all Gauss points of edge
// Success: result should be equal to the edge length
TEST(CompService, integrateScalarFunctionWithGivenValues1D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.5}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::vector<double> f (e.nGP);
    for (int i = 0; i < e.nGP; ++i)
        f[i] = 1.0;

    double res = integrate(e, f);

    EXPECT_FLOAT_EQ(res, e.getLength());
}

// Test of integration of the vector function through 30 deg rotated edge 
// Input data:
// - edge (0;0) - (0.5*sqrt(3), 0.5)
// - test function (returns numvector<double, dimPh>(1.0) in all points) defined as std::function
// Success: result should be equal to the edge length
TEST(CompService, integrateVectorFunction1D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.5}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::function<numvector<double, dimPh>(const Point &)> f = [&](const Point& p){ return numvector<double, dimPh>(1.0); };

    numvector<double, dimPh> res = integrate(e, f);
    numvector<double, dimPh> ref(e.getLength());

    EXPECT_TRUE(equalNumvectors(res, ref));
}

// Test of integration of the vector function through 30 deg rotated edge 
// Input data:
// - edge (0;0) - (0.5*sqrt(3), 0.5)
// - test function (returns numvector<double, dimPh>(1.0) in all points) values in all Gauss points of edge
// Success: result should be equal to the edge length
TEST(CompService, integrateVectorFunctionWithGivenValues1D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({0.866025403785, 0.5}));
    std::vector<std::shared_ptr<Point>> nodes = {n1, n2};

    Edge e (nodes);
    std::vector<numvector<double, dimPh>> f (e.nGP);
    for (int i = 0; i < e.nGP; ++i)
        f[i] = numvector<double, dimPh>(1.0);

    numvector<double, dimPh> res = integrate(e, f);
    numvector<double, dimPh> ref(e.getLength());

    EXPECT_TRUE(equalNumvectors(res, ref));
}

// Test of integration of the scalar function through cell 
// Input data:
// - cell (quad & triag)
// - test function (returns 1.0 in all points) defined as std::function
// Success: result should be equal to the edge length
TEST(CompService, integrateScalarFunction2D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({ 0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({ 0.866025403785, 0.5}));
    std::shared_ptr<Point> n3 = std::make_shared<Point>(Point({ 0.5, 0.866025403785}));
    std::shared_ptr<Point> n4 = std::make_shared<Point>(Point({-0.5, 0.866025403785}));

    std::vector<std::shared_ptr<Point>> nodesQuad  = {n1, n2, n3, n4};
    std::vector<std::shared_ptr<Point>> nodesTriag = {n1, n2, n3};

    std::vector<std::shared_ptr<Edge>> emptyEdges = {};
    Cell cellQuad  (nodesQuad, emptyEdges);
    Cell cellTriag (nodesTriag, emptyEdges);

    std::function<double(const Point &)> f = [&](const Point& p){ return 1.0; };

    double resQuad = integrate(cellQuad, f);
    double resTriag = integrate(cellTriag, f);

    EXPECT_FLOAT_EQ(resQuad,  cellQuad.getArea());
    EXPECT_FLOAT_EQ(resTriag, cellTriag.getArea());
}

// Test of integration of the scalar function through cell 
// Input data:
// - cell (quad & triag)
// - test function (returns 1.0 in all points) defined as std::function
// Success: result should be equal to the edge length
TEST(CompService, integrateVectorFunction2D) {

    std::shared_ptr<Point> n1 = std::make_shared<Point>(Point({ 0.0, 0.0}));
    std::shared_ptr<Point> n2 = std::make_shared<Point>(Point({ 0.866025403785, 0.5}));
    std::shared_ptr<Point> n3 = std::make_shared<Point>(Point({ 0.5, 0.866025403785}));
    std::shared_ptr<Point> n4 = std::make_shared<Point>(Point({-0.5, 0.866025403785}));

    std::vector<std::shared_ptr<Point>> nodesQuad  = {n1, n2, n3, n4};
    std::vector<std::shared_ptr<Point>> nodesTriag = {n1, n2, n3};

    std::vector<std::shared_ptr<Edge>> emptyEdges = {};
    Cell cellQuad  (nodesQuad, emptyEdges);
    Cell cellTriag (nodesTriag, emptyEdges);

    std::function<numvector<double, dimPh>(const Point &)> f = [&](const Point& p){ return numvector<double, dimPh>(1.0); };

    numvector<double, dimPh> resQuad = integrate(cellQuad, f);
    numvector<double, dimPh> resTriag = integrate(cellTriag, f);

    numvector<double, dimPh> refQuad(cellQuad.getArea());
    numvector<double, dimPh> refTriag(cellTriag.getArea());

    EXPECT_TRUE(equalNumvectors(resQuad, refQuad));
    EXPECT_TRUE(equalNumvectors(resTriag, refTriag));
}


}  // namespace