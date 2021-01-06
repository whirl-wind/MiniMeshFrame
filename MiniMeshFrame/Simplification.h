#ifndef SIMPLIFICATION_H_
#define SIMPLIFICATION_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include "DifferentialGeometry.h"
#include <vector>
#include <queue>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialLLT<SparseMatrixType> Solve_LLT;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;

class Simplification
{
public:
	Simplification();
	~Simplification();
	void LoadDDG(DifferentialGeometry* ddg);
	virtual void simplify();
	virtual void simplify_one_step();

	struct edge_v {
		float cost;
		MyMesh::EdgeHandle edge;
		Eigen::Vector4f v;
		edge_v(float c, MyMesh::EdgeHandle e, Eigen::Vector4f vv) {
			cost = c;
			edge = e;
			v = vv;
		}
		friend bool operator < (const edge_v& a, const edge_v& b) {
			return a.cost < b.cost;
		}
		friend bool operator > (const edge_v& a, const edge_v& b) {
			return a.cost > b.cost;
		}
	};
	struct cmp
	{
		bool operator () (const edge_v* a, const edge_v* b) const
		{
			return  a->cost > b->cost;
		}
	};//从小到大

	std::priority_queue<edge_v*, vector<edge_v*>, cmp> q;
protected:
	MyMesh* ptr_mesh_;
	DifferentialGeometry* DDG;
};

class Sp_QEM :
	public Simplification 
{
public:
	Sp_QEM(MyMesh* ptr_mesh_, int Nv, float  Cboundary);
	~Sp_QEM();
	
	void Init();
	bool Qualified();
	float update_min_cost();

	void simplify();
	void simplify_one_step();

	int N_collapse;

private:
	int Nv;
	float Cboundary;
	OpenMesh::FPropHandleT<Eigen::Matrix4f> QEM_f;
	OpenMesh::VPropHandleT<Eigen::Matrix4f> QEM_v;
	OpenMesh::EPropHandleT<Eigen::Matrix4f> QEM_e;
	OpenMesh::EPropHandleT<edge_v*> ev;
};

#endif // !SIMPLIFICATION_H_