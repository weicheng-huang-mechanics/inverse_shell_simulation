#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_thickness, double m_Possion, double m_dt)
{
	YoungM = m_YoungM;
	density = m_density;
	thickness = m_thickness;
	Possion = m_Possion;
	dt = m_dt;

	alpha = YoungM * Possion / ( 1 - Possion * Possion );
	beta = YoungM / ( 2 * (1 + Possion) );

	setupGeometry();

	ndof = 3 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(3 * i + 0) = v_nodes[i](0);
		x(3 * i + 1) = v_nodes[i](1);
		x(3 * i + 2) = v_nodes[i](2);
	}
	x0 = x;

	readInputTriangular();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	for (int i = 0; i < triangularNum; i++)
	{
		int index1 = v_triangularElement[i].nv_1;
		int index2 = v_triangularElement[i].nv_2;
		int index3 = v_triangularElement[i].nv_3;

		massArray(3 * index1 + 0) = massArray(3 * index1 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index1 + 1) = massArray(3 * index1 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index1 + 2) = massArray(3 * index1 + 2) + v_triangularElement[i].area * thickness * density / 3;
		
		massArray(3 * index2 + 0) = massArray(3 * index2 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index2 + 1) = massArray(3 * index2 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index2 + 2) = massArray(3 * index2 + 2) + v_triangularElement[i].area * thickness * density / 3;

		massArray(3 * index3 + 0) = massArray(3 * index3 + 0) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index3 + 1) = massArray(3 * index3 + 1) + v_triangularElement[i].area * thickness * density / 3;
		massArray(3 * index3 + 2) = massArray(3 * index3 + 2) + v_triangularElement[i].area * thickness * density / 3;
	}
}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector3d position, int k)
{
	isConstrained[3 * k + 0] = 1;
	isConstrained[3 * k + 1] = 1;
	isConstrained[3 * k + 2] = 1;

	// Store in the constrained dof vector
	x(3 * k + 0) = position(0);
	x(3 * k + 1) = position(1);
	x(3 * k + 2) = position(2);
}

void elasticPlate::setConstraint(double position, int k)
{
	isConstrained[k] = 1;

	// Store in the constrained dof vector
	x(k) = position;
}

void elasticPlate::setOneVertexBoundaryCondition(double position, int i, int k)
{
	isConstrained[3 * i + k] = 1;

	// Store in the constrained dof vector
	x(3 * i + k) = position;
}

Vector3d elasticPlate::getVertex(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x(3 * i + 0);
	xCurrent(1) = x(3 * i + 1);
	xCurrent(2) = x(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVertexOld(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x0(3 * i + 0);
	xCurrent(1) = x0(3 * i + 1);
	xCurrent(2) = x0(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVelocity(int i)
{
	Vector3d uCurrent;

	uCurrent(0) = ( x(3 * i + 0) - x0(3 * i + 0) ) / dt;
	uCurrent(1) = ( x(3 * i + 1) - x0(3 * i + 1) ) / dt;
	uCurrent(2) = ( x(3 * i + 2) - x0(3 * i + 2) ) / dt;

	uCurrent(0) = u(3 * i + 0);
	uCurrent(1) = u(3 * i + 1);
	uCurrent(2) = u(3 * i + 2);

	return uCurrent;
}

void elasticPlate::readInputTriangular()
{
	triangularNum = 0;
	v_triangularElement.clear();

	for (int i = 0; i < triangular.size(); i++)
	{
		Vector3i triangularCurrent = triangular[i];

		triangularElement m_triangularElement;

		m_triangularElement.nv_1 = triangularCurrent(0);
		m_triangularElement.nv_2 = triangularCurrent(1);
		m_triangularElement.nv_3 = triangularCurrent(2);

		m_triangularElement.x_1 = getVertex(m_triangularElement.nv_1);
		m_triangularElement.x_2 = getVertex(m_triangularElement.nv_2);
		m_triangularElement.x_3 = getVertex(m_triangularElement.nv_3);

		Vector3d e_1 = m_triangularElement.x_2 - m_triangularElement.x_1;
		Vector3d e_2 = m_triangularElement.x_3 - m_triangularElement.x_1;

		m_triangularElement.area = 0.5 * ( e_1.cross(e_2) ).norm();

		m_triangularElement.abar(0, 0) = e_1.dot(e_1);
		m_triangularElement.abar(0, 1) = e_1.dot(e_2);
		m_triangularElement.abar(1, 0) = e_2.dot(e_1);
		m_triangularElement.abar(1, 1) = e_2.dot(e_2);

		m_triangularElement.abarinv = m_triangularElement.abar.inverse();

		// compute b bar
		Vector3i cfaceidx = mesh.F.row(i);
		Vector3d x_temp_1 = getVertex(cfaceidx(0));
		Vector3d x_temp_2 = getVertex(cfaceidx(1));
		Vector3d x_temp_3 = getVertex(cfaceidx(2));

		Matrix3d qc;
		qc.col(0) = x_temp_1;
		qc.col(1) = x_temp_2;
		qc.col(2) = x_temp_3;

		Vector3d edge1 = x_temp_2 - x_temp_1;
		Vector3d edge2 = x_temp_3 - x_temp_1;

		Vector3d cNormal = edge1.cross(edge2);

		Matrix3d oppNormals;
		oppNormals.setZero();

		for (int ii = 0; ii < 3; ii++)
		{
			int edge_idx = mesh.FE(i, ii);
			int edge_orient = mesh.FEorient(i, ii);
			Vector2i edge_vert = mesh.EV.row(edge_idx);

			if (edge_orient == 1)
			{
				int temp = edge_vert(0);
				edge_vert(0) = edge_vert(1);
				edge_vert(1) = temp;
			}

			int oppvert = mesh.Fbendvert(i,3+ii); 

			Vector3d v2 = getVertex( edge_vert(1) );
			Vector3d v3 = getVertex( edge_vert(0) );

			if (oppvert != -1) 
			{
				Vector3d v1 = getVertex(oppvert);

				Vector3d oppNormal_temp = faceNormal(v1, v2, v3);

            	oppNormals.col(ii) = oppNormal_temp;

        	}
        	else if (mesh.cEOpp(edge_idx) != 0) 
        	{
            	Vector3d vec = v3 - v2;
            	Vector3d sysaxi = Vector3d::Zero();

            	int cEOppVal = mesh.cEOpp(edge_idx);

            	if (cEOppVal <= 3) 
            	{
                	sysaxi(cEOppVal - 1) = 1;
            	} 
            	else 
            	{
                	sysaxi(cEOppVal - 4) = -1;
            	}

            	Vector3d oppNormal = sysaxi.cross(vec.normalized()) * 1e5;
            	oppNormals.col(ii) = oppNormal;
        	}
		}

		//cout << " cNormal " << cNormal.transpose() << endl;
		//cout << " oppNormals " << endl;
		//cout << oppNormals << endl;

		Vector3d II = Vector3d::Zero();

    	for (int ii = 0; ii < 3; ii++)
    	{
        	Vector3d oppNormal = oppNormals.col(ii);
        	Vector3d mvec = oppNormal + cNormal;        // averaged normal
        	double mnorms = mvec.norm();                // norm of averaged normal

       	 	int ip1 = (ii + 1) % 3;
        	int ip2 = (ii + 2) % 3;

        	Vector3d qvec = qc.col(ip1) + qc.col(ip2) - 2.0 * qc.col(ii);

        	// Compute II[i]: projection of normal difference onto diagonal vector
        	Vector3d diff_vec = qc.col(ip1) + qc.col(ip2) - 2.0 * qc.col(ii);
        	II(ii) = diff_vec.dot(oppNormal) / mnorms;
        }

        m_triangularElement.bbar(0,0) = II(0)+II(1);
        m_triangularElement.bbar(0,1) = II(0);
        m_triangularElement.bbar(1,0) = II(0);
        m_triangularElement.bbar(1,1) = II(0)+II(2);

        //cout << " bbar " << endl;
        //cout << m_triangularElement.bbar << endl;

		m_triangularElement.arrayIndex = VectorXi::Zero(9);

		m_triangularElement.arrayIndex(0) = 3 * m_triangularElement.nv_1 + 0;
		m_triangularElement.arrayIndex(1) = 3 * m_triangularElement.nv_1 + 1;
		m_triangularElement.arrayIndex(2) = 3 * m_triangularElement.nv_1 + 2;

		m_triangularElement.arrayIndex(3) = 3 * m_triangularElement.nv_2 + 0;
		m_triangularElement.arrayIndex(4) = 3 * m_triangularElement.nv_2 + 1;
		m_triangularElement.arrayIndex(5) = 3 * m_triangularElement.nv_2 + 2;

		m_triangularElement.arrayIndex(6) = 3 * m_triangularElement.nv_3 + 0;
		m_triangularElement.arrayIndex(7) = 3 * m_triangularElement.nv_3 + 1;
		m_triangularElement.arrayIndex(8) = 3 * m_triangularElement.nv_3 + 2;

		v_triangularElement.push_back(m_triangularElement);

		triangularNum = triangularNum + 1;		
	}
}

void elasticPlate::prepareForIteration()
{
	;
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	nv = 0;
	triangularNum = 0;

	ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	v_nodes.clear();
	double a, b, c;
	while(inFile1 >> a >> b >> c)
	{
		Vector3d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		xCurrent(2) = c;

		v_nodes.push_back(xCurrent);
	}
	nv = v_nodes.size();
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/triangleInput.txt");
	v_triangularElement.clear();
	int aa, bb, cc;
	while(inFile2 >> aa >> bb >> cc)
	{
		aa = aa;
		bb = bb;
		cc = cc;

		Vector3i triangularCurrent;
		triangularCurrent(0) = aa;
		triangularCurrent(1) = bb;
		triangularCurrent(2) = cc;
		triangular.push_back(triangularCurrent);
	}
	triangularNum = triangular.size();
	inFile2.close();

    F = MatrixXi::Zero(triangularNum, 3);

    for(int i = 0; i < triangular.size(); i++)
    {
        Vector3i triangularCurrent = triangular[i];

        F.row(i) = triangularCurrent;
    }

    buildMeshConnectivity(F);

    std::cout << "nfaces = " << mesh.nfaces << "\n";
    std::cout << "nedges = " << mesh.nedges << "\n";
    std::cout << "EV:\n" << mesh.EV << "\n";
    std::cout << "EF:\n" << mesh.EF << "\n";
    std::cout << "EOpp:\n" << mesh.EOpp << "\n";
    std::cout << "FE:\n" << mesh.FE << "\n";
    std::cout << "FEorient:\n" << mesh.FEorient << "\n";
    std::cout << "Fbendvert:\n" << mesh.Fbendvert << "\n";
    std::cout << "Fbenddof:\n" << mesh.Fbenddof << "\n";
    std::cout << "cEOpp:\n" << mesh.cEOpp << "\n";

    totalEdge = mesh.nedges;
}

void elasticPlate::buildMeshConnectivity(MatrixXi F) 
{
		mesh.F = F;

    	std::map<std::pair<int, int>, Eigen::Vector2i > edgeFaces;

        int nfaces = (int)F.rows();
        mesh.nfaces = nfaces;

        
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int v0 = F(i, (j + 1) % 3);
                int v1 = F(i, (j + 2) % 3);
                int idx = 0;
                if (v0 > v1)
                {
                    std::swap(v0, v1);
                    idx = 1;
                }

                std::pair<int, int> p(v0, v1);
                auto it = edgeFaces.find(p);
                if (it == edgeFaces.end())
                {
                    edgeFaces[p][idx] = i;
                    edgeFaces[p][1 - idx] = -1;
                }
                else
                {
                    edgeFaces[p][idx] = i;
                }
            }
        }

        

        int nedges = (int)edgeFaces.size();
        mesh.nedges = nedges;

        mesh.FE.resize(nfaces, 3);
        mesh.FEorient.resize(nfaces, 3);
        mesh.EV.resize(nedges, 2);
        mesh.EF.resize(nedges, 2);
        mesh.EOpp.resize(nedges, 2);
        std::map<std::pair<int, int>, int> edgeIndices;

        int idx = 0;
        for (auto it : edgeFaces)
        {
            edgeIndices[it.first] = idx;
            mesh.EV(idx, 0) = it.first.first;
            mesh.EV(idx, 1) = it.first.second;

            if (it.second[0] < it.second[1] )
            {
            	mesh.EF(idx, 0) = it.second[0];
            	mesh.EF(idx, 1) = it.second[1];
            }
            else
            {
            	mesh.EF(idx, 0) = it.second[1];
            	mesh.EF(idx, 1) = it.second[0];
            }

            if (mesh.EF(idx, 0) < 0)
            {
            	int temp = mesh.EF(idx, 0);
            	mesh.EF(idx, 0) = mesh.EF(idx, 1);
            	mesh.EF(idx, 1) = temp;
            }
       
            idx++;
        }

        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int v0 = F(i, (j + 1) % 3);
                int v1 = F(i, (j + 2) % 3);
                if (v0 > v1) std::swap(v0, v1);
                mesh.FE(i, j) = edgeIndices[std::pair<int, int>(v0, v1)];
            }
        }

        for (int i = 0; i < nedges; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                mesh.EOpp(i, j) = oppositeVertex(i, j);
            }
        }

        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int edgeID = mesh.FE(i, j);
        		Vector2i edgeVerts = mesh.EV.row(edgeID);
        
        		int v1 = mesh.F(i, (j+1) % 3);
        		int v2 = mesh.F(i, (j+2) % 3);

        		if ( v1 == edgeVerts(0) && v2 == edgeVerts(1) )
        		{
        			mesh.FEorient(i, j) = 0;
        		}
        		else
        		{
        			mesh.FEorient(i, j) = 1;
        		}
            }
        }

        mesh.Fbendvert = MatrixXi::Constant(nfaces, 6, -1);
    	mesh.Fbenddof  = MatrixXi::Constant(nfaces, 18, -1);

    	for (int i = 0; i < nfaces; ++i) 
    	{
        	std::vector<int> verts = {F(i, 0), F(i, 1), F(i, 2)};

        	for (int j = 0; j < 3; ++j) 
        	{
            	mesh.Fbendvert(i, j) = verts[j];
            	mesh.Fbenddof(i, j * 3 + 0) = 3 * verts[j] + 0;
            	mesh.Fbenddof(i, j * 3 + 1) = 3 * verts[j] + 1;
            	mesh.Fbenddof(i, j * 3 + 2) = 3 * verts[j] + 2;
        	}

        	for (int j = 0; j < 3; ++j) 
        	{
            	int edgeID = mesh.FE(i, j);
            	std::vector<int> eOppVerts = {mesh.EOpp(edgeID, 0), mesh.EOpp(edgeID, 1)};
            	// Remove any that are in face verts
            	for (int v : verts) 
            	{
                	eOppVerts.erase(std::remove(eOppVerts.begin(), eOppVerts.end(), v), eOppVerts.end());
           	 	}
            	if (!eOppVerts.empty()) 
            	{
                	int opp = eOppVerts[0];
                	mesh.Fbendvert(i, 3 + j) = opp;

                	mesh.Fbenddof(i, 9 + j * 3 + 0) = 3 * opp + 0;
               	 	mesh.Fbenddof(i, 9 + j * 3 + 1) = 3 * opp + 1;
                	mesh.Fbenddof(i, 9 + j * 3 + 2) = 3 * opp + 2;
            }
        }
    }

        

    mesh.cEOpp = VectorXi::Zero(mesh.nedges);
}

Vector3d elasticPlate::faceNormal(const Vector3d& q0, const Vector3d& q1, const Vector3d& q2)
{
	Vector3d n = (q1 - q0).cross(q2 - q0);

	return n;
}

int elasticPlate::edgeVertex(int edge, int vertidx) 
{ 
	return mesh.EV(edge, vertidx); 
}

int elasticPlate::edgeFace(int edge, int faceidx)
{ 
	return mesh.EF(edge, faceidx); 
}

int elasticPlate::edgeOppositeVertex(int edge, int faceidx) 
{ 
	return mesh.EOpp(edge, faceidx); 
}

int elasticPlate::faceVertex(int face, int vertidx)
{ 
	return mesh.F(face, vertidx); 
}

int elasticPlate::faceEdge(int face, int vertidx)
{ 
	return mesh.FE(face, vertidx);
} 

int elasticPlate::faceEdgeOrientation(int face, int vertidx)
{
	return mesh.FEorient(face, vertidx);
} 

int elasticPlate::oppositeVertexIndex(int edge, int faceidx)
{
    int face = edgeFace(edge, faceidx);

    if (face == -1)
    {
        return -1;
    }

    for (int j = 0; j < 3; j++)
    {
        if (F(face, j) != edgeVertex(edge, 0) && F(face, j) != edgeVertex(edge, 1))
        return j;
    }

    return -1;
}

int elasticPlate::oppositeVertex(int edge, int faceidx)
{
    int face = edgeFace(edge, faceidx);
    int idx = oppositeVertexIndex(edge, faceidx);

    if (idx == -1)
    {
        return -1;
    }

    return F(face, idx);
}

int elasticPlate::vertexOppositeFaceEdge(int face, int vertidx)
{
    int edge = faceEdge(face, vertidx);
    int edgeorient = faceEdgeOrientation(face, vertidx);

    return edgeOppositeVertex(edge, 1 - edgeorient);
}
