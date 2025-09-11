#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct triangularElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	Vector3d x_1;
	Vector3d x_2;
	Vector3d x_3;

	Matrix2d abar;
	Matrix2d abarinv;

	Matrix2d bbar;

	double area;

	VectorXi arrayIndex;
};


struct Mesh 
{
    int nfaces;
    int nedges;

    MatrixXi F;          // (nfaces x 3)
    MatrixXi EV;         // (nedges x 2) unique edges (vmin, vmax)
    MatrixXi EF;         // (nedges x 2) incident faces, -1 if none
    MatrixXi EOpp;       // (nedges x 2) opposite vertex in each incident face, -1 if none
    MatrixXi FE;         // (nfaces x 3) edge id of each face edge
    MatrixXi FEorient;   // (nfaces x 3) 0 if (vj,vj+1)==EV order, else 1
    MatrixXi Fbendvert;  // (nfaces x 6)
    MatrixXi Fbenddof;   // (nfaces x 18)
    VectorXi cEOpp;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_thickness, double m_Possion, double m_dt);
	~elasticPlate();

	double YoungM;
	double thickness;
	double Possion;
	double dt;
	double density;

	double alpha;
	double beta;

	Vector3d getVertex(int i);
	Vector3d getVertexOld(int i);
	Vector3d getVelocity(int i);

	void setOneVertexBoundaryCondition(double position, int i, int k);

	VectorXd x;
	VectorXd x0;
	VectorXd u;

	std::vector<Vector3d> v_nodes;
	std::vector<triangularElement> v_triangularElement;

	int nv;
	int triangularNum;
	int bendingNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector3d position, int k);
	void setConstraint(double position, int k);

	void readInputEdge();
	void readInputTriangular();

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();

    std::vector<Vector3d> nodes;
    std::vector<Vector3i> triangular;

    void buildMeshConnectivity(MatrixXi F); 
    Mesh mesh;

    Vector3d faceNormal(const Vector3d& q0, const Vector3d& q1, const Vector3d& q2);
    MatrixXi F;

    int totalEdge;

    int oppositeVertexIndex(int edge, int faceidx);
    int oppositeVertex(int edge, int faceidx);
    int vertexOppositeFaceEdge(int face, int vertidx);

    int edgeVertex(int edge, int vertidx); 
    int edgeFace(int edge, int faceidx);
    int edgeOppositeVertex(int edge, int faceidx);

    int faceVertex(int face, int vertidx);
	int faceEdge(int face, int vertidx);
	int faceEdgeOrientation(int face, int vertidx);

	private:
};

#endif
