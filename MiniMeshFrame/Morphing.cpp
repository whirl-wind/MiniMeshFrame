#include "Morphing.h"

Morphing::Morphing()
{
	ptr_mesh_ = new MyMesh();
	ptr_target_mesh_ = new MyMesh();
	DDG = new DifferentialGeometry();
	t = 0.01;
}

Morphing::Morphing(MyMesh* mesh_, MyMesh* target_, float tt)
{
	ptr_mesh_ = mesh_;
	ptr_target_mesh_ = target_;
	t = tt;
}

Morphing::~Morphing()
{
}

void Morphing::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Morphing::morphing(float t_now)
{
}

void Morphing::Init()
{
}

Mp_ARAP::Mp_ARAP(MyMesh* mesh_, MyMesh* target_, DifferentialGeometry* ddg, float tt)
{
	ptr_mesh_ = mesh_;
	ptr_target_mesh_ = target_;
	t = tt;
	LoadDDG(ddg);
	//ptr_backup_ = new MyMesh();
	//ptr_backup_->assign(*ptr_mesh_);
	Init();
}

Mp_ARAP::~Mp_ARAP()
{
}

void Mp_ARAP::morphing(float t_now, int i)
{
	LocalSetL();
	GlobalSolveU(solver, t_now, i);
}

void Mp_ARAP::Init()
{
	int nV = ptr_mesh_->n_vertices();
	int nF = ptr_target_mesh_->n_faces();
	// set two fixed vertices
	if (nV < 2) {
		printf("ERROR::ARAP::Init:\n"
			"\t""need more vertices\n");
		return;
	}
	size_t a = DDG->Boundaries[0].size() / 2;
	auto v1 = ptr_mesh_->from_vertex_handle(DDG->Boundaries[0][0]);

	fixed_vertices_.push_back(v1.idx());
	fixed_position_.push_back(Eigen::Vector3f(ptr_mesh_->point(v1)[0], ptr_mesh_->point(v1)[1], ptr_mesh_->point(v1)[2]));
	L.resize(nF);
	S.resize(nF);

	std::vector<T> A_Triplet;
	SparseMatrixType A(nV, nV);
	for (MyMesh::VertexIter v_it = ptr_target_mesh_->vertices_begin(); v_it != ptr_target_mesh_->vertices_end(); ++v_it) {
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			A_Triplet.push_back(T(i, i, 1));
			continue;
		}
		double sum = 0;

		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_target_mesh_->voh_begin(*v_it); voh_it != ptr_target_mesh_->voh_end(*v_it); ++voh_it) {
			double coe_ij = 0;
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_target_mesh_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			if (!ptr_target_mesh_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_target_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += DDG->cot[t][index];
			}
			// Right Triangle
			triangle = ptr_target_mesh_->face_handle(ptr_target_mesh_->opposite_halfedge_handle(*voh_it));
			if (!ptr_target_mesh_->is_boundary(ptr_target_mesh_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_target_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += DDG->cot[t][index == 0 ? 2 : index - 1];
			}
			// Set
			sum += coe_ij;
			int v_index = ptr_target_mesh_->to_vertex_handle(*voh_it).idx();
			if (v_index != fixed_vertices_[0])
				A_Triplet.push_back(T(i, v_index, -coe_ij));
		}
		A_Triplet.push_back(T(i, i, sum));
	}
	A.setFromTriplets(A_Triplet.begin(), A_Triplet.end());
	solver.compute(A);
}

void Mp_ARAP::LocalSetL()
{
	for (MyMesh::FaceIter f_it = ptr_target_mesh_->faces_begin(); f_it != ptr_target_mesh_->faces_end(); ++f_it)
	{
		int idx = (*f_it).idx();
		Eigen::Vector3d p_(0, 0, 0), q_(0, 0, 0);
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_target_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
			for (int i = 0; i < 3; i++) {
				p_(i) += ptr_mesh_->point(*fv_it)[i];
				q_(i) += ptr_target_mesh_->point(*fv_it)[i];
			}
		}
		int n = x_tmp.size();
		auto p = (ptr_mesh_->point(x_tmp[1]) - ptr_mesh_->point(x_tmp[0])) % (ptr_mesh_->point(x_tmp[2]) - ptr_mesh_->point(x_tmp[0]));
		auto q = (ptr_target_mesh_->point(x_tmp[1]) - ptr_target_mesh_->point(x_tmp[0])) % (ptr_target_mesh_->point(x_tmp[2]) - ptr_target_mesh_->point(x_tmp[0]));
		Eigen::Matrix3d X, U;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < n - 1; j++) {
				X(i, j) = (ptr_mesh_->point(x_tmp[j])[i] - ptr_mesh_->point(x_tmp[j + 1])[i]);
				U(i, j) = (ptr_target_mesh_->point(x_tmp[j])[i] - ptr_target_mesh_->point(x_tmp[j + 1])[i]);
			}
			X(i, n - 1) = p[i];
			U(i, n - 1) = q[i];
		}
		Eigen::Matrix3d J = U * X.inverse();
		//std::cout << J * X - U << std::endl << std::endl;
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
		L[idx] = svd.matrixU() * svd.matrixV().transpose();
		S[idx] = svd.matrixV() * svd.matrixU().transpose() * J;
		//std::cout <<idx<< std::endl << J * X - U << std::endl << std::endl;
	}
}

void Mp_ARAP::GlobalSolveU(Solve_LLT& solver, float t_now, int times)
{
	size_t nV = ptr_mesh_->n_vertices();
	Eigen::VectorXd b_x = Eigen::VectorXd::Zero(nV);
	Eigen::VectorXd b_y = Eigen::VectorXd::Zero(nV);
	Eigen::VectorXd b_z = Eigen::VectorXd::Zero(nV);

	Eigen::Matrix3d I;
	I.setIdentity();

	for (auto v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++)
	{
		int i = (*v_it).idx();
		if (i == fixed_vertices_[0])
		{
			b_x(i) = fixed_position_[0][0];
			b_y(i) = fixed_position_[0][1];
			b_z(i) = fixed_position_[0][2];
			continue;
		}
		/*
		else if(i == fixed_vertices_[1])
		{
			b_x(i) = fixed_position_[1][0];
			b_y(i) = fixed_position_[1][1];
			b_z(i) = fixed_position_[1][2];
			continue;
		}*/
		//auto diff_0 = DDG->ptr_mesh_->point(*v_it);
		Eigen::Vector3d sum = Eigen::Vector3d::Zero();

		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_mesh_->voh_begin(*v_it); voh_it != ptr_mesh_->voh_end(*v_it); ++voh_it)
		{
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_mesh_->face_handle(*voh_it);
			size_t tr = 0;
			int index;
			double coe_ij = 0;
			MyMesh::Point diff_x = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*voh_it));
			MyMesh::Point diff_1(0, 0, 0), diff_2(0, 0, 0);
			if (!ptr_mesh_->is_boundary(*voh_it))
			{
				tr = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				diff_1 = ptr_mesh_->point(*v_it) - diff_x;
				coe_ij += DDG->cot[tr][index];
				sum += DDG->cot[tr][index] * ((1.0 - t_now) * I + t_now * L[tr]) * ((1.0 - t_now) * I + t_now * S[tr]) * Eigen::Vector3d(diff_1[0], diff_1[1], diff_1[2]);
			}
			// Right Triangle
			triangle = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(*voh_it));
			if (!ptr_mesh_->is_boundary(ptr_mesh_->opposite_halfedge_handle(*voh_it)))
			{
				tr = triangle.idx();
				index = 0;
				MyMesh::VertexHandle v2;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				diff_2 = ptr_mesh_->point(*v_it) - diff_x;
				coe_ij += DDG->cot[tr][index == 0 ? 2 : index - 1];
				sum += DDG->cot[tr][index == 0 ? 2 : index - 1] * ((1.0 - t_now) * I + t_now * L[tr]) * ((1.0 - t_now) * I + t_now * S[tr]) * Eigen::Vector3d(diff_2[0], diff_2[1], diff_2[2]);
			}
			//std::cout << diff_1 << std::endl<< diff_2 << std::endl << std::endl;
			// Set
			int j = ptr_mesh_->to_vertex_handle(*voh_it).idx();
			if (j == fixed_vertices_[0])
			{
				b_x(i) -= (-coe_ij) * fixed_position_[0][0];
				b_y(i) -= (-coe_ij) * fixed_position_[0][1];
				b_z(i) -= (-coe_ij) * fixed_position_[0][2];
				continue;
			}
			/*
			else if (i == fixed_vertices_[1])
			{
				b_x(i) -= (-coe_ij) * fixed_position_[1][0];
				b_y(i) -= (-coe_ij) * fixed_position_[1][1];
				b_z(i) -= (-coe_ij) * fixed_position_[1][2];
				continue;
			}*/
		}

		b_x(i) += sum(0);
		b_y(i) += sum(1);
		b_z(i) += sum(2);

	}
	//if (nV <= 10)  std::cout << b_x << std::endl;
	Eigen::VectorXd u_x = solver.solve(b_x);
	Eigen::VectorXd u_y = solver.solve(b_y);
	Eigen::VectorXd u_z = solver.solve(b_z);
	//if (nV <= 10) std::cout << u_x << std::endl;

	MyMesh::Point fix = ptr_target_mesh_->point(ptr_target_mesh_->vertex_handle(fixed_vertices_[0]));
	float dx = (fix[0] - fixed_position_[0][0]) / float(int(1 / t) - 1) * times;
	float dy = (fix[1] - fixed_position_[0][1]) / float(int(1 / t) - 1) * times;
	float dz = (fix[2] - fixed_position_[0][2]) / float(int(1 / t) - 1) * times;

	for (auto v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = MyMesh::Point(u_x(i) + dx, u_y(i) + dy, u_z(i) + dz);
		ptr_mesh_->set_point(*v_it, new_vert);
	}
}
