#include "elasticBendingForce.h"

elasticBendingForce::elasticBendingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
}

elasticBendingForce::~elasticBendingForce()
{
	;
}

void elasticBendingForce::computeFb()
{
    for (int kkk = 0; kkk < plate->triangular.size(); kkk++)
    {
        double thickness = plate->thickness;

        Matrix2d abars = plate->v_triangularElement[kkk].abar;
        Matrix2d bbars = plate->v_triangularElement[kkk].bbar;

        double coeff = pow(thickness, 3) / 12;
        Matrix2d abarinv = plate->v_triangularElement[kkk].abarinv;

        Matrix<double, 4, 18 > bderiv_0;

        Matrix2d b_0 = secondFundamentalForm( kkk, &bderiv_0, 0 );

        Matrix<double, 4, 18 > bderiv;
        Matrix2d b = secondFundamentalForm( kkk, &bderiv, 1 );

        Matrix2d M = abarinv * (b - bbars);
        double dA = 0.5 * sqrt(abars.determinant());

        double alpha, beta;
        alpha = plate->alpha;
        beta = plate->beta;

        double StVK = 0.5 * alpha * pow(M.trace(), 2) + beta * (M * M).trace();
        double result = coeff * dA * StVK;

        Matrix2d temp = alpha * M.trace() * abarinv + 2 * beta * M * abarinv;
        VectorXd derivative;
        derivative = coeff * dA * bderiv_0.transpose() * Map<Vector4d>(temp.data());


        MatrixXd hessian;

        Matrix<double, 1, 18 > inner_0 = bderiv_0.transpose() * Map<Vector4d>(abarinv.data());
        Matrix<double, 1, 18 > inner   = bderiv.transpose() *   Map<Vector4d>(abarinv.data());

        hessian = alpha * inner_0.transpose() * inner;

        Matrix<double, 1, 18 > inner00_0 = abarinv(0, 0) * bderiv_0.row(0) + abarinv(0, 1) * bderiv_0.row(2);
        Matrix<double, 1, 18 > inner01_0 = abarinv(0, 0) * bderiv_0.row(1) + abarinv(0, 1) * bderiv_0.row(3);
        Matrix<double, 1, 18 > inner10_0 = abarinv(1, 0) * bderiv_0.row(0) + abarinv(1, 1) * bderiv_0.row(2);
        Matrix<double, 1, 18 > inner11_0 = abarinv(1, 0) * bderiv_0.row(1) + abarinv(1, 1) * bderiv_0.row(3);
        
        Matrix<double, 1, 18 > inner00 = abarinv(0, 0) * bderiv.row(0) + abarinv(0, 1) * bderiv.row(2);
        Matrix<double, 1, 18 > inner01 = abarinv(0, 0) * bderiv.row(1) + abarinv(0, 1) * bderiv.row(3);
        Matrix<double, 1, 18 > inner10 = abarinv(1, 0) * bderiv.row(0) + abarinv(1, 1) * bderiv.row(2);
        Matrix<double, 1, 18 > inner11 = abarinv(1, 0) * bderiv.row(1) + abarinv(1, 1) * bderiv.row(3);
        

        hessian += 2 * beta *  inner00_0.transpose() * inner00;
        hessian += 2 * beta * (inner01_0.transpose() * inner10 + inner10_0.transpose() * inner01);
        hessian += 2 * beta *  inner11_0.transpose() * inner11;

        hessian = hessian * coeff * dA;

        /*

        Matrix<double, 1, 18 > inner = bderiv.transpose() * Map<Vector4d>(abarinv.data());

        MatrixXd hessian;
        hessian = alpha * inner.transpose() * inner;
        Matrix2d Mainv = M * abarinv;
        for (int i = 0; i < 4; ++i) // iterate over Mainv and abarinv as if they were vectors
        {
            hessian += (alpha * M.trace() * abarinv(i) + 2 * beta * Mainv(i)) * bhess[i];
        }

        Matrix<double, 1, 18 > inner00 = abarinv(0, 0) * bderiv.row(0) + abarinv(0, 1) * bderiv.row(2);
        Matrix<double, 1, 18 > inner01 = abarinv(0, 0) * bderiv.row(1) + abarinv(0, 1) * bderiv.row(3);
        Matrix<double, 1, 18 > inner10 = abarinv(1, 0) * bderiv.row(0) + abarinv(1, 1) * bderiv.row(2);
        Matrix<double, 1, 18 > inner11 = abarinv(1, 0) * bderiv.row(1) + abarinv(1, 1) * bderiv.row(3);
        hessian += 2 * beta * inner00.transpose() * inner00;
        hessian += 2 * beta * (inner01.transpose() * inner10 + inner10.transpose() * inner01);
        hessian += 2 * beta * inner11.transpose() * inner11;

        hessian = hessian * coeff * dA;

        */



        VectorXi arrayNum = plate->mesh.Fbenddof.row(kkk);

        for (int j = 0; j < 18; j++)
        {
            if ( arrayNum(j) > 0 )
            {
                stepper->addForce(arrayNum(j), derivative(j) );
            }

        }

        for (int j = 0; j < 18; j++)
        {
            for (int k = 0; k < 18; k++)
            {
                if ( arrayNum(j) > 0 && arrayNum(k) > 0 )
                {
                    stepper->addJacobian(arrayNum(j), arrayNum(k), hessian(j,k) );
                }
                
            }
        }

    }


}

void elasticBendingForce::computeJb()
{
    ;
}

void elasticBendingForce::setFirstJacobian()
{
    for (int i = 0; i < plate->triangular.size(); i++)
    {
        VectorXi arrayNum = plate->mesh.Fbenddof.row(i);

        for (int j = 0; j < 18; j++)
        {
            for (int k = 0; k < 18; k++)
            {

                if ( arrayNum(j) > 0 && arrayNum(k) > 0 )
                {
                    stepper->addJacobian(arrayNum(j), arrayNum(k), 1);
                }
                
            }
        }
    }
}

Matrix3d elasticBendingForce::crossMat(Vector3d a)
{
    Matrix3d b;

    b<<0,-a(2),a(1),
    a(2),0,-a(0),
    -a(1),a(0),0;

    return b;
}

Vector3d elasticBendingForce::faceNormal(
        const Vector3d& q0,
        const Vector3d& q1,
        const Vector3d& q2,
        Matrix<double, 3, 9>* derivative)
    {
        if (derivative)
        {
            derivative->setZero();
        }

        Vector3d n = (q1 - q0).cross(q2 - q0);

        if (derivative)
        {
            derivative->block(0, 0, 3, 3) += crossMat(q2 - q1);
            derivative->block(0, 3, 3, 3) += crossMat(q0 - q2);
            derivative->block(0, 6, 3, 3) += crossMat(q1 - q0);
        }

        return n;
    }

Matrix2d elasticBendingForce::secondFundamentalForm(
        int kk,
        Matrix<double, 4, 18>* derivative,
        int type)
    {
        if (derivative)
        {
            derivative->resize(4, 18);
            derivative->setZero();
        }


        Matrix<double, 3, 18> IIderiv;
        std::vector < Matrix<double, 18, 18> > IIhess;

        Vector3d II = secondFundamentalFormEntries(kk,  &IIderiv, type);

        Matrix2d result;
        result << II[0] + II[1], II[0], II[0], II[0] + II[2];

        if (derivative)
        {
            derivative->row(0) += IIderiv.row(0);
            derivative->row(0) += IIderiv.row(1);

            derivative->row(1) += IIderiv.row(0);
            derivative->row(2) += IIderiv.row(0);

            derivative->row(3) += IIderiv.row(0);
            derivative->row(3) += IIderiv.row(2);
        }

        return result;
    }

Vector3d elasticBendingForce::secondFundamentalFormEntries(
        int kk,
        Matrix<double, 3, 18>* derivative,
        int type)
    {
        if (derivative)
        {
            derivative->setZero();
        }

        Vector3d II;

        Vector3d oppNormals[3];
        oppNormals[0].setZero();
        oppNormals[1].setZero();
        oppNormals[2].setZero();

        Matrix<double, 3, 9> dn[3];

        dn[0].setZero();
        dn[1].setZero();
        dn[2].setZero();

        Matrix<double, 3, 9> dcn;

        Vector3i cfaceidx = plate->mesh.F.row(kk);
        Vector3d qi, qj, qk;

        if (type == 0)
        {
            qi = plate->getVertex_start(cfaceidx(0));
            qj = plate->getVertex_start(cfaceidx(1));
            qk = plate->getVertex_start(cfaceidx(2));
        }


        if (type == 1)
        {
            qi = plate->getVertex(cfaceidx(0));
            qj = plate->getVertex(cfaceidx(1));
            qk = plate->getVertex(cfaceidx(2));
        }
        

        Matrix3d qc;
        qc.col(0) = qi;
        qc.col(1) = qj;
        qc.col(2) = qk;

        Vector3d cNormal = faceNormal(qi, qj, qk, &dcn );

        for (int ii = 0; ii < 3; ii++)
        {
            int edge_idx = plate->mesh.FE(kk, ii);
            int edge_orient = plate->mesh.FEorient(kk, ii);
            Vector2i edge_vert = plate->mesh.EV.row(edge_idx);

            if (edge_orient == 1)
            {
                int temp = edge_vert(0);
                edge_vert(0) = edge_vert(1);
                edge_vert(1) = temp;
            }

            int oppvert = plate->mesh.Fbendvert(kk,3+ii); 

            Vector3d v2, v3;

            if (type == 0)
            {
                v2 = plate->getVertex_start( edge_vert(1) );
                v3 = plate->getVertex_start( edge_vert(0) );
            }

            if (type == 1)
            {
                v2 = plate->getVertex( edge_vert(1) );
                v3 = plate->getVertex( edge_vert(0) );
            }

            if (oppvert != -1) 
            {
                Vector3d v1;

                if (type == 0)
                {
                    v1 = plate->getVertex_start(oppvert);
                }

                if (type == 1)
                {
                    v1 = plate->getVertex(oppvert);
                }

                oppNormals[ii] = faceNormal(v1,v2,v3,  &dn[ii]);
            }

            else if (plate->mesh.cEOpp(edge_idx) != 0) 
            {
                Vector3d vec = v3 - v2;
                Vector3d sysaxi = Vector3d::Zero();

                int cEOppVal = plate->mesh.cEOpp(edge_idx);

                if (cEOppVal <= 3) 
                {
                    sysaxi(cEOppVal - 1) = 1;
                } 
                else 
                {
                    sysaxi(cEOppVal - 4) = -1;
                }

                Vector3d oppNormal = sysaxi.cross(vec.normalized()) * 1e5;
                oppNormals[ii] = oppNormal;
            }
        }

        Vector3d qs[3];
        Vector3d mvec[3];
        Vector3d qvec[3];
        double mnorms[3];

        qs[0] = qi;
        qs[1] = qj;
        qs[2] = qk;

        for (int i = 0; i < 3; i++)
        {
            mvec[i] = oppNormals[i] + cNormal;
            mnorms[i] = mvec[i].norm();
        }


        for (int i = 0; i < 3; i++)
        {
            int ip1 = (i + 1) % 3;
            int ip2 = (i + 2) % 3;
            qvec[i] = qs[ip1] + qs[ip2] - 2.0 * qs[i];
        }

        for (int i = 0; i < 3; i++)
        {
            int ip1 = (i + 1) % 3;
            int ip2 = (i + 2) % 3;
            II[i] = (qs[ip1] + qs[ip2] - 2.0 * qs[i]).dot(oppNormals[i]) / mnorms[i];
            if (derivative)
            {
                derivative->block<1, 3>(i, 3 * i) += -2.0 * oppNormals[i].transpose() / mnorms[i];
                derivative->block<1, 3>(i, 3 * ip1) += 1.0 * oppNormals[i].transpose() / mnorms[i];
                derivative->block<1, 3>(i, 3 * ip2) += 1.0 * oppNormals[i].transpose() / mnorms[i];

                derivative->block<1, 3>(i, 9 + 3 * i) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 0);
                derivative->block<1, 3>(i, 3 * ip2) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 3);
                derivative->block<1, 3>(i, 3 * ip1) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 6);

                derivative->block<1, 3>(i, 9 + 3 * i) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 0); // 9+3i是i的对顶点，也就是qi_bar
                derivative->block<1, 3>(i, 3 * ip2) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 3);   // ip2是i在三角内部左侧的点，表明dn中第二个索引应该组装到qk
                derivative->block<1, 3>(i, 3 * ip1) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 6);

                derivative->block<1, 3>(i, 0) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 0);
                derivative->block<1, 3>(i, 3) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 3);
                derivative->block<1, 3>(i, 6) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 6);
            }
            
        }

        return II;
    }
