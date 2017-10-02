#include "SbplSplineMotionPrimitives.hpp"
#include <iostream> //FIXME remove after debug
#include <set>

using namespace base::geometry;
namespace motion_planning_libraries 
{

SbplSplineMotionPrimitives::SbplSplineMotionPrimitives() {}
  
SbplSplineMotionPrimitives::SbplSplineMotionPrimitives(const SplinePrimitivesConfig& config) :
    config(config)
{
    validateConfig(config);
    
    radPerDiscreteAngle = ((config.assume_bidirectional ? 1.0 : 2.0) *M_PI) / config.numAngles;
    primitivesByAngle.resize(config.numAngles);
    generatePrimitives(config);
}

std::vector<Eigen::Vector2i> SbplSplineMotionPrimitives::generateDestinationCells(const SplinePrimitivesConfig& config) const
{
    const int R = config.destinationCircleRadius;
    const int R2 = R*R;
    std::cout << "generateDestinationCells for R = " << R << ":";
    std::vector<Eigen::Vector2i> result;
    for(int x=-R; x<=R; ++x)
    {
        for(int y=-R; y<=R; ++y)
        {
            if(x*x + y*y <= R*R && (x!=0 || y!=0))
            {
                result.emplace_back(x,y);
                std::cout << " [" << result.back().transpose() << ']';
            }
        }
    }
    std::cout << '\n';
    return result;

    //for cellSkipFactor < 1 the algorithm generates the same coordinates several times
    //therefore we put them in a set :-) This is simpler than writing a more complex algorithm
    std::set<std::pair<int, int>> coordinates; //use pair instead of Vector2i because operator< is implemented for pair
    double r = 0;
    for(int n = 0; r < config.destinationCircleRadius; ++n, r += config.cellSkipFactor)
    {
        //see https://en.wikipedia.org/wiki/Midpoint_circle_algorithm
        if(r <= 0) //happens for cellSkipFactor < 1
            r = 1; 
        int x = r;
        int y = 0;
        int err = 0;

        std::cout << "R: " << r << std::endl;
        
        while (x >= y)
        {
            coordinates.emplace(x,y);
            coordinates.emplace(y,x);
            coordinates.emplace(-y, x);
            coordinates.emplace(-x, y);
            coordinates.emplace(-x, -y);
            coordinates.emplace(-y, -x);
            coordinates.emplace(y, -x);
            coordinates.emplace(x,-y);

            y += 1;
            err += 1 + 2*y;
            if (2*(err-x) + 1 > 0)
            {
                x -= 1;
                err += 1 - 2*x;
            }
        }
    }
    return result;
}

void SbplSplineMotionPrimitives::generatePrimitives(const SplinePrimitivesConfig& config)
{
    std::vector<Eigen::Vector2i> destinationCells = generateDestinationCells(config);
  
    for(int startAngle = 0; startAngle < config.numAngles; ++startAngle)
    {
        generatePrimitivesForAngle(startAngle, destinationCells);
    }
    
    //debug code below
    size_t count = 0;
    for(int angle = 0; angle < config.numAngles; ++angle)
    {
        count += getPrimitiveForAngle(angle).size();
    }
}

void SbplSplineMotionPrimitives::generatePrimitivesForAngle(const int startAngle,
                                                            std::vector<Eigen::Vector2i> destinationCells)
{
    /* Main idea:
     * For each destination cell: generate primitives from (0,0) to that cell.
     * The amount of primitives is defined by the config.
     * Additionally generate point turn primitives from (0,0) to (0,0)  
     */
    
    //NOTE since the destinationCells form a circle, we do not need to rotate them.
    int id = 0; //the "same" primitives should have the same id for each start angle according to sbpl
    const double radStartAngle = startAngle * radPerDiscreteAngle;
    const double epsilon = 0.1; //makes the distinction between forward/backward/lateral easier
    const std::vector<int> endAngles = generateEndAngles(startAngle, config);
        
    for(const Eigen::Vector2i& dest : destinationCells)
    { 
        //rotate destination for easier checks. If it is rotated we just need to check
        //whether destRot.x is above 0 (epsilon)  to know if this is a forward or backward motion
        const Eigen::Vector2d destRot = Eigen::Rotation2D<double>(-radStartAngle) * dest.cast<double>();

        //forward and backward movements
        for(int endAngle : endAngles)
        {
            //forward movement
            if(config.generateForwardMotions &&
               destRot.x() > epsilon)
            {
                SplinePrimitive prim = getPrimitive(startAngle, endAngle, dest,
                                                    id, SplinePrimitive::SPLINE_MOVE_FORWARD);
                primitivesByAngle[startAngle].push_back(prim);
                ++id;
            }
            //backward movement
            else if(config.generateBackwardMotions &&
                    destRot.x() < - epsilon)
            {
                SplinePrimitive prim = getPrimitive(startAngle, endAngle,
                                                    dest, id, SplinePrimitive::SPLINE_MOVE_BACKWARD);
                prim.startAngleRad = startAngle * radPerDiscreteAngle;
                prim.endAngleRad = endAngle * radPerDiscreteAngle;
                primitivesByAngle[startAngle].push_back(prim);
                ++id;
            }
        }
        
        //lateral movement
        if(config.generateLateralMotions &&
           destRot.x() <=  epsilon &&
           destRot.x() >=  -epsilon)
        {
            //the robot is driving sidewards. spline is calculated as if the robot is rotated 90° 
            // and moving forward. angles are fixed afterwards
            int rotatedStartAngle = 0;
            if(destRot.y() < 0)
            {//right rotate
                rotatedStartAngle = (startAngle - config.numAngles / 4) % config.numAngles;
            }
            else
            {
                rotatedStartAngle = (startAngle + config.numAngles / 4) % config.numAngles;
            }
            
            SplinePrimitive prim = getPrimitive(rotatedStartAngle, rotatedStartAngle,
                                                dest, id, SplinePrimitive::SPLINE_MOVE_LATERAL);
            prim.startAngle = startAngle;
            prim.endAngle = startAngle;
            prim.endAngleRad = radStartAngle;
            prim.startAngleRad = radStartAngle;
            primitivesByAngle[startAngle].push_back(prim);
            ++id;
        }
    }
    
    if(config.generatePointTurnMotions)
    {
        //point turns are a special case
        for(int angle = 0; angle < config.numAngles; ++ angle)
        {
            if(angle == startAngle)
                continue;
            
            SplinePrimitive prim;
            prim.startAngle = startAngle;
            prim.endAngle = angle;
            prim.startAngleRad = radStartAngle;
            prim.endAngleRad = angle * radPerDiscreteAngle;
            prim.id = id;
            prim.motionType = SplinePrimitive::SPLINE_POINT_TURN;
            prim.endPosition << 0, 0;
            primitivesByAngle[startAngle].push_back(prim);
            ++id;
        }
    }
    
}

SplinePrimitive SbplSplineMotionPrimitives::getPrimitive(const int startAngle,
                                                         const int endAngle,
                                                         const Eigen::Vector2i destination,
                                                         const int primId,
                                                         const SplinePrimitive::Type& type) const
{
    SplinePrimitive prim;
    const bool backward = (type == SplinePrimitive::SPLINE_MOVE_BACKWARD);
    const double radStartAngle = startAngle * radPerDiscreteAngle + (backward ? M_PI : 0.0);
    const base::Vector2d start(0, 0); 
    const base::Vector2d startDirection = start + Eigen::Rotation2D<double>(radStartAngle) * Eigen::Vector2d::UnitX();
    
    // if we move bidirectional, make sure that actual endAngle has less than 90 degrees difference to start angle
    const int actual_endAngle = config.assume_bidirectional ? (startAngle + (endAngle - startAngle + config.numAngles/2) % config.numAngles - config.numAngles/2) : endAngle;

    const double radEndAngle = actual_endAngle * radPerDiscreteAngle + (backward ? M_PI : 0.0);
    const base::Vector2d end(destination.cast<double>() * config.gridSize);
    const base::Vector2d endDirection = end + Eigen::Rotation2D<double>(radEndAngle) * Eigen::Vector2d::UnitX();
    
    std::vector<base::Vector2d> points{start, startDirection, end, endDirection};
    std::vector<SplineBase::CoordinateType> types{SplineBase::ORDINARY_POINT,
                                                    SplineBase::TANGENT_POINT_FOR_PRIOR,
                                                    SplineBase::ORDINARY_POINT,
                                                    SplineBase::TANGENT_POINT_FOR_PRIOR};
                                                    
    prim.spline = Spline2(config.splineGeometricResolution, config.splineOrder);
    prim.spline.interpolate(points, std::vector<double>(), types);
    prim.startAngle = startAngle;
    prim.endAngleRad = radEndAngle;
    prim.startAngleRad = radStartAngle;
    prim.endAngle = endAngle;
    prim.endPosition = destination;
    prim.id = primId;
    prim.motionType = type;

    return prim;
}

std::vector<int> SbplSplineMotionPrimitives::generateEndAngles(const int startAngle, const SplinePrimitivesConfig& config) const
{
    //otherwise the calculation below gets more complicated
    assert(config.assume_bidirectional || config.numAngles % 2 == 0);
//    const int actual_angles = (config.assume_bidirectional ? 2 : 1) * config.numAngles;
    std::cout << "End angles for start angle = " << startAngle << ": ";
    
    //the left/right distinction is necessary in case of an even number of end angles (FIXME which actually is supposed to be forbidden?)
    const int numAnglesLeftSide = (config.numEndAngles - 1) / 2;
    const int numAnglesRightSide = (config.numEndAngles - (config.numAngles % 2)) / 2;
    std::vector<int> result;

    for(int angle = startAngle - numAnglesRightSide; angle <= startAngle + numAnglesLeftSide; ++angle)
    {
        result.push_back((angle+config.numAngles) % config.numAngles);
        std::cout << result.back() << ", ";
    }

    std::cout << '\n';
    return result;
}


const std::vector<SplinePrimitive>& SbplSplineMotionPrimitives::getPrimitiveForAngle(const int angle) const
{
    assert(angle >= 0);
    assert(angle < config.numAngles);
    return primitivesByAngle[angle];
}

const SplinePrimitivesConfig& SbplSplineMotionPrimitives::getConfig() const
{
    return config;
}

void SbplSplineMotionPrimitives::validateConfig(const SplinePrimitivesConfig& config) const
{
    if(config.gridSize <= 0)
        throw std::runtime_error("gridSize has to be > 0");
    
    if(config.numAngles <= 0)
        throw std::runtime_error("numAngles has to be > 0");
    
    if(config.numEndAngles <= 0)
        throw std::runtime_error("numEndAngles has to be > 0");
    
    if(config.destinationCircleRadius <= 0)
        throw std::runtime_error("destinationCircleRadius has to be > 0");
    
    if(config.numAngles % 2 != 0)
        throw std::runtime_error("numAngles has to be even");
        
    if(config.numEndAngles > config.numAngles / (config.assume_bidirectional ? 1 : 2))
        throw std::runtime_error("numEndAngles has to be <= numAngles / 2");
    
    if(config.splineGeometricResolution <= 0)
        throw std::runtime_error("splineGeometricResolution has to be > 0");
    
    if(config.splineOrder < 3)
        throw std::runtime_error("splineOrder has to be >= 3");
}



}//end namespace motion_planning_libraries
