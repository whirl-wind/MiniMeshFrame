#include "Smoothing.h"
#include <stdlib.h>
#include <time.h> 
#include <cmath>

Smoothing::Smoothing()
{
	ptr_mesh_ = nullptr;
}

Smoothing::~Smoothing()
{
}

void Smoothing::AddNoise(double scale)
{
	if (scale == 0 || ptr_mesh_ == nullptr) return;
	if (ptr_mesh_->n_vertices() == 0) return;

	int randmax = 5000;
	srand((unsigned)time(NULL));

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		int r = rand() % randmax;
		double x = r / (randmax - 1.0) * scale;
		r = rand() % randmax;
		double y = r / (randmax - 1.0) * scale;
		r = rand() % randmax;
		double z = r / (randmax - 1.0) * scale;

		//if (ptr_mesh_->is_boundary(*v_it)) continue;
		
		MyMesh::Point new_point = ptr_mesh_->point(*v_it) + MyMesh::Point(x,y,z);
		ptr_mesh_->set_point(*v_it, new_point);
	}
	ptr_mesh_->update_normals();
}

void Smoothing::LoadMesh(MyMesh* ptr_)
{
	ptr_mesh_ = ptr_;
}

void Smoothing::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Smoothing::SmoothMesh()
{
}

BilateralMeshDenoising::BilateralMeshDenoising()
{
	ptr_mesh_ = nullptr;
	iter_num = 5;
	is_global_ = true;
	sig_s = 0.5;
	sig_c = 0.5;
	n_lam = 0.5;
	p_lam = 0.5;
}

BilateralMeshDenoising::~BilateralMeshDenoising()
{
}

void BilateralMeshDenoising::Set_LocalParameter(int num, double sig)
{
	iter_num = num;
	sig_s = sig;
}

void BilateralMeshDenoising::Set_GlobalParameter(double n_l, double p_l, double sig)
{
	n_lam = n_l;
	p_lam = p_l;
	sig_s = sig;
}

void BilateralMeshDenoising::SmoothMesh()
{
	if (!ptr_mesh_) return;
	if (ptr_mesh_->n_faces() < 2) return;
	if (!DDG) {
		DDG = new DifferentialGeometry();
		DDG->update_mesh(ptr_mesh_);
	}
	double V0 = DDG->Caculate_Volume();
	if(is_global_) Global_EnergyMinimization();
	else Local_Iterative();
	double V1 = DDG->Caculate_Volume();
	ScaleMesh(sqrt(V0 / V1));
	ptr_mesh_->update_normals();
}

void BilateralMeshDenoising::ScaleMesh(double rate)
{
	MyMesh::Point o = MyMesh::Point(DDG->Geometry_Analysis_Center(0), DDG->Geometry_Analysis_Center(1), DDG->Geometry_Analysis_Center(2));
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		MyMesh::Point p = ptr_mesh_->point((*v_it));
		MyMesh::Point v = rate * p + (1 - rate) * o;
		ptr_mesh_->set_point(*v_it, v);
	}
}

void BilateralMeshDenoising::SetGlobal(bool global)
{
	is_global_ = global;
}

void BilateralMeshDenoising::update_NormalField()
{
	fnormal.resize(ptr_mesh_->n_faces());
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		int i = (*f_it).idx();
		fnormal[i].setZero();
		Eigen::Vector3d ni = Eigen::Vector3d(ptr_mesh_->normal(f_it)[0], ptr_mesh_->normal(f_it)[1], ptr_mesh_->normal(f_it)[2]);
		sig_c = 0;
		int count = 0;
		for (auto ff_it = ptr_mesh_->ff_begin(f_it.handle()); ff_it != ptr_mesh_->ff_end(f_it.handle()); ++ff_it) {
			int j = (*ff_it).idx();
			sig_c += (DDG->FaceCenter[i] - DDG->FaceCenter[j]).norm();
			count++;
		}
		sig_c = sig_c / (1.0 * count);
		for (auto ff_it = ptr_mesh_->ff_begin(f_it.handle()); ff_it != ptr_mesh_->ff_end(f_it.handle()); ++ff_it) {
			int j = (*ff_it).idx();
			Eigen::Vector3d nj = Eigen::Vector3d(ptr_mesh_->normal(ff_it)[0], ptr_mesh_->normal(ff_it)[1], ptr_mesh_->normal(ff_it)[2]);
			double Ws = exp(-(ni - nj).squaredNorm() / (2 * sig_s * sig_s));
			double Wc = exp(-(DDG->FaceCenter[i] - DDG->FaceCenter[j]).squaredNorm() / (2 * sig_c * sig_c));
			fnormal[i] += DDG->FaceArea[j] * Ws * Wc * nj;
		}
		fnormal[i].normalize();
	}
}

void BilateralMeshDenoising::update_VertexPosition()
{
	if (is_global_) {
		int row, col;
		row = ptr_mesh_->n_faces() * 3 * 3;
		col = ptr_mesh_->n_vertices() * 3;

		Eigen::VectorXd b = Eigen::VectorXd(col);
		for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
			int i = (*v_it).idx();
			b(i * 3 + 0) = ptr_mesh_->point(*v_it)[0];
			b(i * 3 + 1) = ptr_mesh_->point(*v_it)[1];
			b(i * 3 + 2) = ptr_mesh_->point(*v_it)[2];
		}
		SparseMatrixType D;
		SparseMatrixType N;
		vector<T> tripletList_D;
		vector<T> tripletList_N;
		D.resize(row, col);
		N.resize(row / 3, row);

		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
		{
			int loop = 0;
			int k = (*f_it).idx();

			tripletList_N.push_back(T(k + 0, (k + 0) * 3 + 0, fnormal[k](0)));
			tripletList_N.push_back(T(k + 0, (k + 0) * 3 + 1, fnormal[k](1)));
			tripletList_N.push_back(T(k + 0, (k + 0) * 3 + 2, fnormal[k](2)));
			tripletList_N.push_back(T(k + 1, (k + 1) * 3 + 0, fnormal[k](0)));
			tripletList_N.push_back(T(k + 1, (k + 1) * 3 + 1, fnormal[k](1)));
			tripletList_N.push_back(T(k + 1, (k + 1) * 3 + 2, fnormal[k](2)));
			tripletList_N.push_back(T(k + 2, (k + 2) * 3 + 0, fnormal[k](0)));
			tripletList_N.push_back(T(k + 2, (k + 2) * 3 + 1, fnormal[k](1)));
			tripletList_N.push_back(T(k + 2, (k + 2) * 3 + 2, fnormal[k](2)));

			auto fv_it = ptr_mesh_->fv_begin(f_it.handle());
			int  begin = (*fv_it).idx();
			int i = begin;
			int j;
			++fv_it;
			double f_los;
			for (; fv_it != ptr_mesh_->fv_end(f_it.handle()); ++fv_it) {
				j = (*fv_it).idx();
				tripletList_D.push_back(T((k + loop) * 3 + 0, i * 3 + 0,  1));
				tripletList_D.push_back(T((k + loop) * 3 + 1, i * 3 + 1,  1));
				tripletList_D.push_back(T((k + loop) * 3 + 2, i * 3 + 2,  1));
				tripletList_D.push_back(T((k + loop) * 3 + 0, j * 3 + 0, -1));
				tripletList_D.push_back(T((k + loop) * 3 + 1, j * 3 + 1, -1));
				tripletList_D.push_back(T((k + loop) * 3 + 2, j * 3 + 2, -1));
				i = j;
				loop++;
			}
			j = begin;
			tripletList_D.push_back(T((k + loop) * 3 + 0, i * 3 + 0, 1));
			tripletList_D.push_back(T((k + loop) * 3 + 1, i * 3 + 1, 1));
			tripletList_D.push_back(T((k + loop) * 3 + 2, i * 3 + 2, 1));
			tripletList_D.push_back(T((k + loop) * 3 + 0, j * 3 + 0, -1));
			tripletList_D.push_back(T((k + loop) * 3 + 1, j * 3 + 1, -1));
			tripletList_D.push_back(T((k + loop) * 3 + 2, j * 3 + 2, -1));
		}
		D.setFromTriplets(tripletList_D.begin(), tripletList_D.end()); //根据三元组列表生成稀疏矩阵
		D.makeCompressed();
		N.setFromTriplets(tripletList_N.begin(), tripletList_N.end()); //根据三元组列表生成稀疏矩阵
		N.makeCompressed();

		SparseMatrixType G = D.transpose() * N.transpose() * N * D;
		SparseMatrixType I;
		I.resize(col, col);
		I.setIdentity();
		
		G = (1 - p_lam) * G + p_lam * I;
		b = p_lam * b;

		Solve_LLT solver;
		solver.compute(G);
		Eigen::VectorXd x = Eigen::VectorXd(col);
		x = solver.solve(b);
	
		for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
			int i = (*v_it).idx();
			ptr_mesh_->set_point(*v_it, MyMesh::Point(x(i * 3 + 0), x(i * 3 + 1), x(i * 3 + 2)));
		}
	}
	else {
		for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
			int i = (*v_it).idx();
			MyMesh::Point vi = ptr_mesh_->point(v_it);
			Eigen::Vector3d xi = Eigen::Vector3d(vi[0], vi[1], vi[2]);
			Eigen::Vector3d di;
			di.setZero();
			int count = 0;
			for (auto vf_it = ptr_mesh_->vf_begin(v_it.handle()); vf_it != ptr_mesh_->vf_end(v_it.handle()); vf_it++) {
				int j = (*vf_it).idx();
				di += (fnormal[j].transpose() * (DDG->FaceCenter[j] - xi))(0) * fnormal[j];
				count++;
			}
			di /= (1.0 * count);
			vi += MyMesh::Point(di(0), di(1), di(2));
			ptr_mesh_->set_point(*v_it, vi);
		}
	}
	ptr_mesh_->update_normals();
}

void BilateralMeshDenoising::Local_Iterative()
{
	for (int i = 0; i < iter_num; i++) {
		DDG->Caculate_FaceArea();
		DDG->Build_FaceCenter();

		update_NormalField();
		update_VertexPosition();

		DDG->Clear_FaceCenter();
	}
}

void BilateralMeshDenoising::Global_EnergyMinimization()
{
	DDG->Caculate_FaceArea();
	DDG->Build_FaceCenter();

	int row, col;
	row = col = ptr_mesh_->n_faces() * 3;
	SparseMatrixType BL;
	SparseMatrixType A;
	vector<T> tripletList;
	vector<T> tripletList_A;
	BL.resize(row, col);
	A.resize(row, col);
	Eigen::VectorXd n0 = Eigen::VectorXd(row);
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		int i = (*f_it).idx();
		double sq_ai = sqrt(DDG->FaceArea[i]);
		tripletList.push_back(T(i * 3 + 0, i * 3 + 0, sq_ai));
		tripletList.push_back(T(i * 3 + 1, i * 3 + 1, sq_ai));
		tripletList.push_back(T(i * 3 + 2, i * 3 + 2, sq_ai));
		tripletList_A.push_back(T(i * 3 + 0, i * 3 + 0, DDG->FaceArea[i]));
		tripletList_A.push_back(T(i * 3 + 1, i * 3 + 1, DDG->FaceArea[i]));
		tripletList_A.push_back(T(i * 3 + 2, i * 3 + 2, DDG->FaceArea[i]));
		Eigen::Vector3d ni = Eigen::Vector3d(ptr_mesh_->normal(f_it)[0], ptr_mesh_->normal(f_it)[1], ptr_mesh_->normal(f_it)[2]);
		
		n0(i * 3 + 0) = ni(0);
		n0(i * 3 + 1) = ni(1);
		n0(i * 3 + 2) = ni(2);

		sig_c = 0;
		int count = 0;
		for (auto ff_it = ptr_mesh_->ff_begin(f_it.handle()); ff_it != ptr_mesh_->ff_end(f_it.handle()); ++ff_it) {
			int j = (*ff_it).idx();
			sig_c += (DDG->FaceCenter[i] - DDG->FaceCenter[j]).norm();
			count++;
		}
		sig_c = sig_c / (1.0 * count);
		for (auto ff_it = ptr_mesh_->ff_begin(f_it.handle()); ff_it != ptr_mesh_->ff_end(f_it.handle()); ++ff_it) {
			int j = (*ff_it).idx();
			Eigen::Vector3d nj = Eigen::Vector3d(ptr_mesh_->normal(ff_it)[0], ptr_mesh_->normal(ff_it)[1], ptr_mesh_->normal(ff_it)[2]);
			double Ws = exp(-(ni - nj).squaredNorm() / (2 * sig_s * sig_s));
			double Wc = exp(-(DDG->FaceCenter[i] - DDG->FaceCenter[j]).squaredNorm() / (2 * sig_c * sig_c));
			double Wij = DDG->FaceArea[j] * Ws * Wc;
			tripletList.push_back(T(i * 3 + 0, j * 3 + 0, -sq_ai * Wij));
			tripletList.push_back(T(i * 3 + 1, j * 3 + 1, -sq_ai * Wij));
			tripletList.push_back(T(i * 3 + 2, j * 3 + 2, -sq_ai * Wij));
		}
	}
	BL.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	BL.makeCompressed();
	A.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
	A.makeCompressed();

	SparseMatrixType G;
	G.resize(row, col);
	G = (1 - n_lam) * BL.transpose() * BL + n_lam * A; // 这里没*2，是因为下面b也没*2


	Solve_LLT solver;
	solver.compute(G);
	Eigen::VectorXd n = Eigen::VectorXd(row);
	Eigen::VectorXd b = Eigen::VectorXd(row);
	b = n_lam * A * n0;
	n = solver.solve(b);
	
	fnormal.resize(ptr_mesh_->n_faces());
	for (int i = 0; i < ptr_mesh_->n_faces(); i++) {
		fnormal[i] = Eigen::Vector3d(n(i * 3 + 0), n(i * 3 + 1), n(i * 3 + 2));
	}

	update_VertexPosition();

	DDG->Clear_FaceCenter();
}
