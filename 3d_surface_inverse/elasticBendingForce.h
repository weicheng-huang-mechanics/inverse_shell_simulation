#ifndef ELASTICBENDINGFORCE_H
#define ELASTICBENDINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticBendingForce
{
public:
	elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticBendingForce();

	void computeFb();
	void computeJb();
    void setFirstJacobian();
    
    VectorXd TotalForceVec;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    Vector3d faceNormal(
            const Vector3d& qi0,
            const Vector3d& qi1,
            const Vector3d& qi2,
            Matrix<double, 3, 9>* derivative);

    Matrix2d secondFundamentalForm(
            int kk,
            Matrix<double, 4, 18>* derivative,
            int type);

    Vector3d secondFundamentalFormEntries(
            int kk,
            Matrix<double, 3, 18>* derivative,
            int type);

    Matrix3d crossMat(Vector3d a);
};

#endif
