#ifndef ELASTICSTRETCHINGFORCE_H
#define ELASTICSTRETCHINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticStretchingForce
{
public:
	elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticStretchingForce();
	void computeFs();
	void computeJs();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    Matrix3d I3;
    Matrix3d Z3;

    void place3(MatrixXd& M, int bi, int bj, Matrix3d B);

};

#endif
