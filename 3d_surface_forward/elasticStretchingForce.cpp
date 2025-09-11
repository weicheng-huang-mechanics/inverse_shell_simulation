#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	I3 << 1,0,0,
          0,1,0,
          0,0,1;

	Z3 << 0,0,0,
          0,0,0,
          0,0,0;
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
	{
		int index1 = plate->v_triangularElement[i].nv_1;
		int index2 = plate->v_triangularElement[i].nv_2;
		int index3 = plate->v_triangularElement[i].nv_3;

		Vector3d vi = plate->getVertex(index1);
		Vector3d vj = plate->getVertex(index2);
		Vector3d vk = plate->getVertex(index3);

		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		Vector3d e_1 = vj - vi;  
		Vector3d e_2 = vk - vi; 

		Matrix2d A;

		A(0,0) = e_1.dot(e_1);
		A(0,1) = e_1.dot(e_2);
		A(1,0) = e_2.dot(e_1);
		A(1,1) = e_2.dot(e_2);

		VectorXd gradA1, gradA2, gradA3, gradA4;
		MatrixXd hessA1, hessA2, hessA3, hessA4;

    	gradA1 = VectorXd::Zero(9);
    	gradA2 = VectorXd::Zero(9);
    	gradA3 = VectorXd::Zero(9);
    	gradA4 = VectorXd::Zero(9);

    	hessA1 = MatrixXd::Zero(9, 9);
    	hessA2 = MatrixXd::Zero(9, 9);
    	hessA3 = MatrixXd::Zero(9, 9);
    	hessA4 = MatrixXd::Zero(9, 9);


    	gradA1.segment(0,3) = - 2 * e_1; 
		gradA1.segment(3,3) =   2 * e_1; 

		gradA2.segment(0,3) = - e_1 - e_2;
		gradA2.segment(3,3) =         e_2;        
		gradA2.segment(6,3) =   e_1;


		gradA3.segment(0,3) = - e_2 - e_1;
		gradA3.segment(3,3) =   e_2;
		gradA3.segment(6,3) =         e_1;

		gradA4.segment(0,3) = - 2 * e_2;
		gradA4.segment(6,3) =   2 * e_2;

	
		place3(hessA1,0,0,  2.0*I3); place3(hessA1,0,1, -2.0*I3); place3(hessA1,0,2, Z3);
		place3(hessA1,1,0, -2.0*I3); place3(hessA1,1,1,  2.0*I3); place3(hessA1,1,2, Z3);
		place3(hessA1,2,0,     Z3 ); place3(hessA1,2,1,     Z3 ); place3(hessA1,2,2, Z3);

		place3(hessA2,0,0,  2.0*I3); place3(hessA2,0,1, -1.0*I3); place3(hessA2,0,2, -1.0*I3);
    	place3(hessA2,1,0, -1.0*I3); place3(hessA2,1,1,     Z3 ); place3(hessA2,1,2,  1.0*I3);
    	place3(hessA2,2,0, -1.0*I3); place3(hessA2,2,1,  1.0*I3); place3(hessA2,2,2,     Z3 );

    	hessA3 = hessA2;

    	place3(hessA4,0,0,  2.0*I3); place3(hessA4,0,1,   Z3    ); place3(hessA4,0,2, -2.0*I3);
    	place3(hessA4,1,0,    Z3  ); place3(hessA4,1,1,   Z3    ); place3(hessA4,1,2,    Z3  );
    	place3(hessA4,2,0, -2.0*I3); place3(hessA4,2,1,   Z3    ); place3(hessA4,2,2,  2.0*I3);


    	double thickness = plate->thickness;
    	double coefficient = thickness / 4;

    	Matrix2d abar = plate->v_triangularElement[i].abar;
    	Matrix2d abarInv = plate->v_triangularElement[i].abarinv;

    	double dA = 0.5 * sqrt( abar.determinant() );

    	double alpha, beta;
    	alpha = plate->alpha;
    	beta = plate->beta;

    	Matrix2d M = abarInv * (A - abar);

    	double StVK = 0.5 * alpha * ( M.trace() ) * ( M.trace() ) + beta * (M * M).trace();


	    Matrix2d tempMatrix;
    	tempMatrix = alpha * M.trace() * abarInv + 2 * beta * M * abarInv;

    	VectorXd derivative;

    	derivative = VectorXd::Zero(9);

    	derivative = dA * coefficient * (gradA1 * tempMatrix(0,0) + gradA2 * tempMatrix(0,1) + gradA3 * tempMatrix(1,0) + gradA4 * tempMatrix(1,1) );  

    	VectorXd inner;

    	inner = (gradA1 * abarInv(0,0) + gradA2 * abarInv(0,1) + gradA3 * abarInv(1,0) + gradA4 * abarInv(1,1) );

    	MatrixXd hessian;

    	hessian = MatrixXd::Zero(9, 9);

    	hessian = (coefficient * dA) * alpha * (inner * inner.transpose());

    	Matrix2d MainV = M * abarInv;

    	double coeff_i;
    
    	coeff_i = alpha * M.trace() * abarInv(0, 0) + 2 * beta * MainV(0, 0);
    	hessian = hessian + ( coefficient * dA ) * coeff_i * hessA1;

    	coeff_i = alpha * M.trace() * abarInv(0, 1) + 2 * beta * MainV(0, 1);
    	hessian = hessian + ( coefficient * dA ) * coeff_i * hessA2;

    	coeff_i = alpha * M.trace() * abarInv(1, 0) + 2 * beta * MainV(1, 0);
    	hessian = hessian + ( coefficient * dA ) * coeff_i * hessA3;

    	coeff_i = alpha * M.trace() * abarInv(1, 1) + 2 * beta * MainV(1, 1);
    	hessian = hessian + ( coefficient * dA ) * coeff_i * hessA4;

    	VectorXd inner001 = abarInv(0,0) * gradA1 + abarInv(0,1) * gradA3; 
		VectorXd inner011 = abarInv(0,0) * gradA2 + abarInv(0,1) * gradA4;
		VectorXd inner101 = abarInv(1,0) * gradA1 + abarInv(1,1) * gradA3;
		VectorXd inner111 = abarInv(1,0) * gradA2 + abarInv(1,1) * gradA4;

		hessian = hessian + (coefficient * dA) * 2 * beta * (inner001 * inner001.transpose());
		hessian = hessian + (coefficient * dA) * 2 * beta * (inner011 * inner101.transpose() + inner101 * inner011.transpose());
		hessian = hessian + (coefficient * dA) * 2 * beta * (inner111 * inner111.transpose());

		for (int j = 0; j < 9; j++)
        {
            stepper->addForce(arrayIndex(j), derivative(j) );
        }


		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), hessian(j, k) );
            }
        }

	}
}

void elasticStretchingForce::computeJs()
{
	;
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
	{
		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), 1);
            }
        }
	}
}

void elasticStretchingForce::place3(MatrixXd& M, int bi, int bj, Matrix3d B) 
{
  M.block<3,3>(3*bi, 3*bj) = B;
}