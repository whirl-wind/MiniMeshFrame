#pragma once
#include <vector>
#include <cmath>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>

typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseLU<SparseMatrixType> Solve_LU;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

using std::vector;
constexpr double MY_PI = 3.1415926;

enum CellType
{
	NoneCell = 0,
	BarycentricCell = 1,
	VoronoiCell = 2,
	MixedVoronoiCell = 3,
};

enum LaplacianType
{
	NoneLaplacian = 0,
	UniformLaplacian = 1,
	CotangentLaplacian = 2,
};

class DifferentialGeometry
{
public:
	DifferentialGeometry();
	~DifferentialGeometry();

	void	update_mesh(MyMesh *ptr);
	void	SetAnalysisCenter(Eigen::Vector3d c);
	void	SetBoundaries();

	void	SetCurrentLaplacian(int i);
	int		GetCurrentLaplacian(int i);
	void	SetCurrentCell(int i);
	int		GetCurrentCell(int i);
	void	Build_UniformLaplacian();
	void	Build_CotangentLaplacian();

	void	Build_BarycentricCell();
	void	Build_VoronoiCell();
	void	Build_MixedVoronoiCell();

	void	Build_FaceCenter();
	void	Clear_FaceCenter();
	void	Caculate_FaceArea();

	void	Caculate_Cot();
	double	Caculate_Volume();
	void	Caculate_MeanCurvature();
	void	Caculate_GaussianCurvature();
	void	Get_MeanCurvature(vector<double> &x);
	void	Get_GaussianCurvature(vector<double> &x);

	std::tuple<double, double, double> BarycentricCoordinate(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d P);

	bool is_build_barycentriccell;
	bool is_build_voronoicell;
	bool is_build_mixedvoronoicell;

	bool is_caculated_meancurvature;
	bool is_caculated_gaussiancurvature;

	double bbox_x_min;
	double bbox_x_max;
	double bbox_y_min;
	double bbox_y_max;
	double bbox_z_min;
	double bbox_z_max;

	double													average_edgelength;
	vector<double>											FaceArea;
	std::vector<std::vector<double>>						cot;
	vector<Eigen::Vector3d>									FaceCenter;
	Eigen::Vector3d											Geometry_Analysis_Center;
	vector< vector<MyMesh::HalfedgeHandle> > 				Boundaries;

	MyMesh* ptr_mesh_;

	vector<double> 											BarycentricCellArea;
	vector<double> 											VoronoiCellArea;
	vector<double> 											MixedVoronoiCellArea;
private:

	CellType					current_cell_type_;
	LaplacianType				current_laplacian_type_;
	
	SparseMatrixType			Uniform_Laplacian; //Lfi = Sum_j (fi -fj)
	SparseMatrixType			Cotangent_Laplacian;
	
	vector<double> 											mean_curvature;
	vector<double> 											gaussian_curvature;
};

