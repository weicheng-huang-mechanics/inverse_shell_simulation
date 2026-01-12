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

        Matrix<double, 4, 18 > bderiv;
        std::vector<Matrix<double, 18 , 18 > > bhess;

        Matrix2d b = secondFundamentalForm( kkk, &bderiv, &bhess );

        Matrix2d M = abarinv * (b - bbars);
        double dA = 0.5 * sqrt(abars.determinant());

        double alpha, beta;
        alpha = plate->alpha;
        beta = plate->beta;

        double StVK = 0.5 * alpha * pow(M.trace(), 2) + beta * (M * M).trace();
        double result = coeff * dA * StVK;

        Matrix2d temp = alpha * M.trace() * abarinv + 2 * beta * M * abarinv;
        VectorXd derivative;
        derivative = coeff * dA * bderiv.transpose() * Map<Vector4d>(temp.data());

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



        VectorXi arrayNum = plate->mesh.Fbenddof.row(kkk);

        for (int j = 0; j < 18; j++)
        {
            if ( arrayNum(j) >= 0 )
            {
                stepper->addForce(arrayNum(j), derivative(j) );
            }

        }

        for (int j = 0; j < 18; j++)
        {
            for (int k = 0; k < 18; k++)
            {
                if ( arrayNum(j) >= 0 && arrayNum(k) >= 0 )
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

                if ( arrayNum(j) >= 0 && arrayNum(k) >= 0 )
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
        Matrix<double, 3, 9>* derivative,
        std::vector<Matrix<double, 9, 9> >* hessian)
    {
        if (derivative)
        {
            derivative->setZero();
        }

        if (hessian)
        {
            hessian->resize(3);
            for (int i = 0; i < 3; i++) 
            {
                (*hessian)[i].setZero();
            }
        }

        Vector3d n = (q1 - q0).cross(q2 - q0);

        if (derivative)
        {
            derivative->block(0, 0, 3, 3) += crossMat(q2 - q1);
            derivative->block(0, 3, 3, 3) += crossMat(q0 - q2);
            derivative->block(0, 6, 3, 3) += crossMat(q1 - q0);
        }

        if (hessian)
        {
            for (int j = 0; j < 3; j++)
            {
                Vector3d ej(0, 0, 0);
                ej[j] = 1.0;
                Matrix3d ejc = crossMat(ej);
                (*hessian)[j].block(0, 3, 3, 3) -= ejc;
                (*hessian)[j].block(0, 6, 3, 3) += ejc;
                (*hessian)[j].block(3, 6, 3, 3) -= ejc;
                (*hessian)[j].block(3, 0, 3, 3) += ejc;
                (*hessian)[j].block(6, 0, 3, 3) -= ejc;
                (*hessian)[j].block(6, 3, 3, 3) += ejc;
            }
        }

        return n;
    }

Matrix2d elasticBendingForce::secondFundamentalForm(
        int kk,
        Matrix<double, 4, 18>* derivative,
        std::vector<Matrix<double, 18, 18> >* hessian)
    {
        if (derivative)
        {
            derivative->resize(4, 18);
            derivative->setZero();
        }

        if (hessian)
        {
            hessian->resize(4);
            for (int i = 0; i < 4; i++)
            {
                (*hessian)[i].resize(18, 18);
                (*hessian)[i].setZero();
            }
        }

        Matrix<double, 3, 18> IIderiv;
        std::vector < Matrix<double, 18, 18> > IIhess;

        Vector3d II = secondFundamentalFormEntries(kk,  &IIderiv , &IIhess );

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
        if (hessian)
        {
            (*hessian)[0] += IIhess[0];
            (*hessian)[0] += IIhess[1];

            (*hessian)[1] += IIhess[0];
            (*hessian)[2] += IIhess[0];

            (*hessian)[3] += IIhess[0];
            (*hessian)[3] += IIhess[2];
        }

        return result;
    }

Vector3d elasticBendingForce::secondFundamentalFormEntries(
        int kk,
        Matrix<double, 3, 18>* derivative,
        std::vector<Matrix<double, 18, 18> >* hessian)
    {
        if (derivative)
        {
            derivative->setZero();
        }

        if (hessian)
        {
            hessian->resize(3);
            for (int i = 0; i < 3; i++)
            {
                (*hessian)[i].setZero();
            }
        }

        Vector3d II;

        Vector3d oppNormals[3];
        oppNormals[0].setZero();
        oppNormals[1].setZero();
        oppNormals[2].setZero();

        Matrix<double, 3, 9> dn[3];
        std::vector<Matrix<double, 9, 9> > hn[3];

        dn[0].setZero();
        dn[1].setZero();
        dn[2].setZero();

        
        for (int i = 0; i < 3; ++i) 
        {
            hn[i].resize(3);

            for (int j = 0; j < 3; ++j) 
            {
                hn[i][j] = Matrix<double, 9, 9>::Zero(); 
            }
        }

        Matrix<double, 3, 9> dcn;
        std::vector<Matrix<double, 9, 9> > hcn;

        Vector3i cfaceidx = plate->mesh.F.row(kk);
        Vector3d qi = plate->getVertex(cfaceidx(0));
        Vector3d qj = plate->getVertex(cfaceidx(1));
        Vector3d qk = plate->getVertex(cfaceidx(2));

        Matrix3d qc;
        qc.col(0) = qi;
        qc.col(1) = qj;
        qc.col(2) = qk;

        Vector3d cNormal = faceNormal(qi, qj, qk,   &dcn , &hcn );

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

            Vector3d v2 = plate->getVertex( edge_vert(1) );
            Vector3d v3 = plate->getVertex( edge_vert(0) );

            if (oppvert != -1) 
            {
                Vector3d v1 = plate->getVertex(oppvert);

                oppNormals[ii] = faceNormal(v1,v2,v3,  &dn[ii] , &hn[ii]);
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
            if (hessian)
            {
                int ip1 = (i + 1) % 3;
                int ip2 = (i + 2) % 3;

                int miidx[3];
                miidx[0] = 9 + 3 * i;
                miidx[1] = 3 * ip2;
                miidx[2] = 3 * ip1;

                Matrix3d dnij[3];
                for (int j = 0; j < 3; j++)
                    dnij[j] = dn[i].block<3, 3>(0, 3 * j);

                for (int j = 0; j < 3; j++)
                {

                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip1) += (1.0 / mnorms[i]) * dnij[j].transpose();
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip2) += (1.0 / mnorms[i]) * dnij[j].transpose();
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * i) += (-2.0 / mnorms[i]) * dnij[j].transpose();

                    Vector3d dnijTm = dnij[j].transpose() * mvec[i];
                    Vector3d dcnjTm = (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]);
                    Matrix3d term3 = dnijTm * oppNormals[i].transpose();

                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;

                    Matrix3d term4 = dcnjTm * oppNormals[i].transpose();

                    (*hessian)[i].block<3, 3>(3 * j, 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
                    (*hessian)[i].block<3, 3>(3 * j, 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
                    (*hessian)[i].block<3, 3>(3 * j, 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;

                    (*hessian)[i].block<3, 3>(3 * ip1, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
                    (*hessian)[i].block<3, 3>(3 * ip2, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
                    (*hessian)[i].block<3, 3>(3 * i, miidx[j]) += (-2.0 / mnorms[i]) * dnij[j];

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dnij[j].transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
                    }

                    Matrix3d term1 = oppNormals[i] * dnijTm.transpose();
                    (*hessian)[i].block<3, 3>(3 * ip1, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
                    (*hessian)[i].block<3, 3>(3 * ip2, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
                    (*hessian)[i].block<3, 3>(3 * i, miidx[j]) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;

                    Matrix3d term2 = oppNormals[i] * dcnjTm.transpose();
                    (*hessian)[i].block<3, 3>(3 * ip1, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
                    (*hessian)[i].block<3, 3>(3 * ip2, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
                    (*hessian)[i].block<3, 3>(3 * i, 3 * j) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;

                    Vector3d dnijTq = dnij[j].transpose() * qvec[i];

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
                    }

                    double qdoto = qvec[i].dot(oppNormals[i]);

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dnij[k];
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dcn.block<3, 3>(0, 3 * k);
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dnij[k];
                        (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dcn.block<3, 3>(0, 3 * k);

                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));

                        for (int l = 0; l < 3; l++)
                        {
                            (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
                            (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hcn[l].block<3, 3>(3 * j, 3 * k);
                            (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (1.0 / mnorms[i]) * qvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
                        }
                    }
                }
            }
        }

        return II;
    }
