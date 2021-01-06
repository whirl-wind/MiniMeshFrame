#include "Deformation.h"

Deformation::Deformation()
{
	ptr_mesh_ = new MyMesh();
	ptr_mesh_->add_property(fixed);
}

Deformation::Deformation(MyMesh* ptr_object_, vector<Interaction*> IXs, Interaction* select_ix)
{
	ptr_mesh_ = ptr_object_;
	ptr_mesh_->add_property(fixed);
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		ptr_mesh_->property(fixed, *v_it) = false;
	}
	for (vector<Interaction*>::iterator ix = IXs.begin();  ix != IXs.end(); ++ix) {
		if ((*ix)->fix) {
			for (vector<MyMesh::VertexHandle>::iterator v_it = (*ix)->point_array.begin(); v_it != (*ix)->point_array.end(); ++v_it) {
				ptr_mesh_->property(fixed, *v_it) = true;
			}
		}
	}
	for (vector<MyMesh::VertexHandle>::iterator v_it = select_ix->point_array.begin(); v_it != select_ix->point_array.end(); ++v_it) {
		ptr_mesh_->property(fixed, *v_it) = true;
	}
}

Deformation::~Deformation()
{
	ptr_mesh_->remove_property(fixed);
}

void Deformation::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Deformation::deform()
{
}

void Deformation::Init()
{
}

void Deformation::update_IXs(vector<Interaction*> IXs, Interaction* select_ix)
{
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		ptr_mesh_->property(fixed, *v_it) = false;
	}
	for (vector<Interaction*>::iterator ix = IXs.begin(); ix != IXs.end(); ++ix) {
		if ((*ix)->fix) {
			for (vector<MyMesh::VertexHandle>::iterator v_it = (*ix)->point_array.begin(); v_it != (*ix)->point_array.end(); ++v_it) {
				ptr_mesh_->property(fixed, *v_it) = true;
			}
		}
	}
	if (select_ix) {
		for (vector<MyMesh::VertexHandle>::iterator v_it = select_ix->point_array.begin(); v_it != select_ix->point_array.end(); ++v_it) {
			ptr_mesh_->property(fixed, *v_it) = true;
		}
	}
}

Df_ARAP::Df_ARAP(MyMesh* ptr_object_, vector<Interaction*> IXs, Interaction* select_ix, DifferentialGeometry* ddg)
{
	ptr_mesh_ = ptr_object_;
	ptr_mesh_->add_property(fixed);
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		ptr_mesh_->property(fixed, *v_it) = false;
	}
	update_IXs(IXs, select_ix);
	LoadDDG(ddg);
	setTimes(5);
	L.resize(ptr_mesh_->n_faces());
}

Df_ARAP::~Df_ARAP()
{
}

void Df_ARAP::deform()
{
	for (int count = 0; count < times; count++)
	{
		LocalSetL();
		GlobalSolveU(solver);
	}
}

void Df_ARAP::Init()
{
	size_t nV = ptr_mesh_->n_vertices();
	SparseMatrixType A(nV, nV);
	GlobalMatrixA(A);
	//if(nV<=10) std::cout << A << std::endl;

	solver.compute(A);

	if (times < 0 || times > 100)
		times = 10;
	//std::cout << "time = " << times << std::endl;
}

void Df_ARAP::update_IXs(vector<Interaction*> IXs, Interaction* select_ix)
{
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		ptr_mesh_->property(fixed, *v_it) = false;
	}
	for (vector<Interaction*>::iterator ix = IXs.begin(); ix != IXs.end(); ++ix) {
		if ((*ix)->fix) {
			for (vector<MyMesh::VertexHandle>::iterator v_it = (*ix)->point_array.begin(); v_it != (*ix)->point_array.end(); ++v_it) {
				ptr_mesh_->property(fixed, *v_it) = true;
			}
		}
	}
	if (select_ix) {
		for (vector<MyMesh::VertexHandle>::iterator v_it = select_ix->point_array.begin(); v_it != select_ix->point_array.end(); ++v_it) {
			ptr_mesh_->property(fixed, *v_it) = true;
		}
	}
}

void Df_ARAP::setTimes(int time)
{
	times = time;
}

void Df_ARAP::LocalSetL()
{
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		int idx = (*f_it).idx();
		Eigen::Vector3d p_(0, 0, 0), q_(0, 0, 0);
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
			for (int i = 0; i < 3; i++) {
				p_(i) += DDG->ptr_mesh_->point(*fv_it)[i];
				q_(i) += ptr_mesh_->point(*fv_it)[i];
			}
		}
		int n = x_tmp.size();
		auto p = (DDG->ptr_mesh_->point(x_tmp[1])- DDG->ptr_mesh_->point(x_tmp[0])) % (DDG->ptr_mesh_->point(x_tmp[2]) - DDG->ptr_mesh_->point(x_tmp[0]));
		auto q = (ptr_mesh_->point(x_tmp[1]) - ptr_mesh_->point(x_tmp[0])) % (ptr_mesh_->point(x_tmp[2]) - ptr_mesh_->point(x_tmp[0]));
		Eigen::Matrix3d X, U;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < n-1; j++) {
				X(i, j) = (DDG->ptr_mesh_->point(x_tmp[j])[i] - DDG->ptr_mesh_->point(x_tmp[j + 1])[i]);
				U(i, j) = (ptr_mesh_->point(x_tmp[j])[i] - ptr_mesh_->point(x_tmp[j + 1])[i]);
			}
			X(i, n - 1) = p[i];
			U(i, n - 1) = q[i];
		}
		Eigen::Matrix3d J = U * X.inverse();
		//std::cout << J * X - U << std::endl << std::endl;
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
		L[idx] = svd.matrixU() * svd.matrixV().transpose();
		//std::cout <<idx<< std::endl << J * X - U << std::endl << std::endl;
	}
}

void Df_ARAP::GlobalSolveU(Solve_LLT& solver)
{
	size_t nV = ptr_mesh_->n_vertices();
	Eigen::VectorXd b_x = Eigen::VectorXd::Zero(nV);
	Eigen::VectorXd b_y = Eigen::VectorXd::Zero(nV);
	Eigen::VectorXd b_z = Eigen::VectorXd::Zero(nV);

	for (auto v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++)
	{
		int i = (*v_it).idx();
		if (ptr_mesh_->property(fixed, *v_it))
		{
			b_x(i) = ptr_mesh_->point(*v_it)[0];
			b_y(i) = ptr_mesh_->point(*v_it)[1];
			b_z(i) = ptr_mesh_->point(*v_it)[2];
			continue;
		}
		auto diff_0 = DDG->ptr_mesh_->point(*v_it);
		Eigen::Vector3d sum = Eigen::Vector3d::Zero();
		
		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_mesh_->voh_begin(*v_it); voh_it != ptr_mesh_->voh_end(*v_it); ++voh_it)
		{
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_mesh_->face_handle(*voh_it);
			size_t tr = 0;
			int index;
			double coe_ij = 0; 
			MyMesh::Point diff_x = DDG->ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*voh_it));
			MyMesh::Point diff_1(0, 0, 0), diff_2(0, 0, 0);
			if (!ptr_mesh_->is_boundary(*voh_it))
			{
				tr = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				diff_1 = DDG->ptr_mesh_->point(*v_it)- diff_x;
				coe_ij += DDG->cot[tr][index];
				sum += DDG->cot[tr][index] * L[tr] * Eigen::Vector3d(diff_1[0], diff_1[1], diff_1[2]);
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
				diff_2 = DDG->ptr_mesh_->point(*v_it) - diff_x;
				coe_ij += DDG->cot[tr][index == 0 ? 2 : index - 1];
				sum += DDG->cot[tr][index == 0 ? 2 : index - 1] * L[tr] * Eigen::Vector3d(diff_2[0], diff_2[1], diff_2[2]);
			}
			//std::cout << diff_1 << std::endl<< diff_2 << std::endl << std::endl;
			// Set
			if (ptr_mesh_->property(fixed, ptr_mesh_->to_vertex_handle(*voh_it)))
			{
				b_x(i) -= (-coe_ij) * ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*voh_it))[0];
				b_y(i) -= (-coe_ij) * ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*voh_it))[1];
				b_z(i) -= (-coe_ij) * ptr_mesh_->point(ptr_mesh_->to_vertex_handle(*voh_it))[2];
			}
			/*else if (j == fixed_vertices_[1])
			{
				b_x(i) -= (-coe_ij) * fixed_coords_[1][0];
				b_y(i) -= (-coe_ij) * fixed_coords_[1][1];
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

	for (auto v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		int i = (*v_it).idx();
		MyMesh::Point new_vert = MyMesh::Point(u_x(i), u_y(i), u_z(i));
		ptr_mesh_->set_point(*v_it, new_vert);
	}
}

void Df_ARAP::GlobalMatrixA(SparseMatrixType& A)
{
	std::vector<T> A_Triplet;
	int nV = ptr_mesh_->n_vertices();

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		int i = (*v_it).idx();
		if (ptr_mesh_->property(fixed, *v_it))
		{
			A_Triplet.push_back(T(i, i, 1));
			continue;
		}
		double sum = 0;

		for (MyMesh::VertexOHalfedgeIter voh_it = ptr_mesh_->voh_begin(*v_it); voh_it != ptr_mesh_->voh_end(*v_it); ++voh_it) {
			double coe_ij = 0;
			// Left Triangle
			MyMesh::FaceHandle triangle = ptr_mesh_->face_handle(*voh_it);
			size_t t = 0;
			int index;
			if (!ptr_mesh_->is_boundary(*voh_it))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += DDG->cot[t][index];
			}
			// Right Triangle
			triangle = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(*voh_it));
			if (!ptr_mesh_->is_boundary(ptr_mesh_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == i) break;
					index++;
				}
				coe_ij += DDG->cot[t][index == 0 ? 2 : index - 1];
			}
			// Set
			sum += coe_ij;
			int v_index = ptr_mesh_->to_vertex_handle(*voh_it).idx();
			if (!ptr_mesh_->property(fixed, ptr_mesh_->to_vertex_handle(*voh_it)))
				A_Triplet.push_back(T(i, v_index, -coe_ij));
		}
		A_Triplet.push_back(T(i, i, sum));
	}
	A.setFromTriplets(A_Triplet.begin(), A_Triplet.end());
}
