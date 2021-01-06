#include "DifferentialGeometry.h"
#include <iostream>

DifferentialGeometry::DifferentialGeometry() {
	ptr_mesh_ = new MyMesh();
	current_cell_type_ = NoneCell;
	current_laplacian_type_ = NoneLaplacian;

	is_build_barycentriccell = false;
	is_build_voronoicell = false;
	is_build_mixedvoronoicell = false;
	is_caculated_meancurvature = false;
	is_caculated_gaussiancurvature = false;
}

DifferentialGeometry::~DifferentialGeometry()
{
	BarycentricCellArea.clear();
	VoronoiCellArea.clear();
	MixedVoronoiCellArea.clear();
	mean_curvature.clear();
	gaussian_curvature.clear();
	FaceArea.clear();
	for (int i = 0; i < cot.size(); i++) {
		cot[i].clear();
	}
	cot.clear();
}

void DifferentialGeometry::update_mesh(MyMesh * ptr)
{
	ptr_mesh_ = ptr;
	current_cell_type_ = NoneCell;
	current_laplacian_type_ = NoneLaplacian;

	is_build_barycentriccell = false;
	is_build_voronoicell = false;
	is_build_mixedvoronoicell = false;
	is_caculated_meancurvature = false;
	is_caculated_gaussiancurvature = false;

	BarycentricCellArea.clear();
	VoronoiCellArea.clear();
	MixedVoronoiCellArea.clear();
	mean_curvature.clear();
	gaussian_curvature.clear();
	FaceArea.clear();

	if (ptr_mesh_->n_faces() == 0) return;

	Caculate_FaceArea();
	Caculate_Cot();

	Build_UniformLaplacian();
	Build_CotangentLaplacian();

	Geometry_Analysis_Center(0) = (bbox_x_max + bbox_x_min) / 2.0;
	Geometry_Analysis_Center(1) = (bbox_y_max + bbox_y_min) / 2.0;
	Geometry_Analysis_Center(2) = (bbox_z_max + bbox_z_min) / 2.0;
	//SetAnalysisCenter(Eigen::Vector3d(((bbox_x_max + bbox_x_min) / 2.0, (bbox_y_max + bbox_y_min) / 2.0, (bbox_z_max + bbox_z_min) / 2.0)));
	SetBoundaries();
}

void DifferentialGeometry::SetAnalysisCenter(Eigen::Vector3d c)
{
	Geometry_Analysis_Center = c;
}

void DifferentialGeometry::SetBoundaries()
{
	Boundaries.clear();
	vector<bool> check_in_boundary;
	check_in_boundary.resize(ptr_mesh_->n_halfedges());
	for (int i = 0; i < ptr_mesh_->n_halfedges(); i++) {
		check_in_boundary[i] = true;
	}
	for (MyMesh::HalfedgeIter he_it = ptr_mesh_->halfedges_begin(); he_it != ptr_mesh_->halfedges_end(); ++he_it) {
		//glNormal3fv(ptr_mesh_->normal(ptr_mesh_->from_vertex_handle(*he_it)).data());
		if (ptr_mesh_->is_boundary(*he_it)) {
			if (check_in_boundary[(*he_it).idx()]) {
				vector<MyMesh::HalfedgeHandle> new_boundary;
				new_boundary.push_back((*he_it));
				check_in_boundary[(*he_it).idx()] = false;
				MyMesh::HalfedgeHandle he_next = ptr_mesh_->next_halfedge_handle((*he_it));
				while (he_next != (*he_it)) {
					new_boundary.push_back(he_next);
					check_in_boundary[he_next.idx()] = false;
					he_next = ptr_mesh_->next_halfedge_handle(he_next);
				}
				Boundaries.push_back(new_boundary);
			}
		}
	}
}

void DifferentialGeometry::SetCurrentLaplacian(int i)
{
	current_laplacian_type_ = LaplacianType(i);
}

int DifferentialGeometry::GetCurrentLaplacian(int i)
{
	
	return current_laplacian_type_;
}

void DifferentialGeometry::SetCurrentCell(int i)
{
	current_cell_type_ = CellType(i);
	switch (i) {
	case BarycentricCell:
		if (is_build_barycentriccell) return;
		else Build_BarycentricCell();
		break;
	case VoronoiCell:
		if (is_build_voronoicell) return;
		else Build_VoronoiCell();
		break;
	case MixedVoronoiCell:
		if (is_build_mixedvoronoicell) return;
		else Build_MixedVoronoiCell();
		break;
	default: 
		break;
	}
}

int DifferentialGeometry::GetCurrentCell(int i)
{
	return current_cell_type_;
}

void DifferentialGeometry::Build_UniformLaplacian()
{
	int row, col;
	row = col = ptr_mesh_->n_vertices();
	current_laplacian_type_ = UniformLaplacian;

	vector<T> tripletList;
	Uniform_Laplacian.resize(row, col);

	bbox_x_min = bbox_y_min = bbox_z_min = INT_MAX;
	bbox_x_max = bbox_y_max = bbox_z_max = INT_MIN;

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		int idx_1 = v_it.handle().idx();
		double x = ptr_mesh_->point(v_it)[0];
		double y = ptr_mesh_->point(v_it)[1];
		double z = ptr_mesh_->point(v_it)[2];
		if (x < bbox_x_min) bbox_x_min = x;
		else if (x > bbox_x_max) bbox_x_max = x;
		if (y < bbox_y_min) bbox_y_min = y;
		else if (y > bbox_y_max) bbox_y_max = y;
		if (z < bbox_z_min) bbox_z_min = z;
		else if (z > bbox_z_max) bbox_z_max = z;
		
		int count = 0;
		for (MyMesh::VertexVertexIter vv_it = ptr_mesh_->vv_begin(v_it); vv_it != ptr_mesh_->vv_end(v_it); vv_it++) {
			int idx_2 = vv_it.handle().idx();
			tripletList.push_back(T(idx_1, idx_2, 1));
			count++;
		}
		tripletList.push_back(T(idx_1, idx_1, -count));
	}

	Uniform_Laplacian.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	Uniform_Laplacian.makeCompressed();

}

void DifferentialGeometry::Build_CotangentLaplacian()
{
	Build_MixedVoronoiCell();
	
	int row, col;
	row = col = ptr_mesh_->n_vertices();
	current_laplacian_type_ = CotangentLaplacian;

	vector<T> tripletList;
	Cotangent_Laplacian.resize(row, col);

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		int idx_1 = v_it.handle().idx();
		double area_2 = 2 * MixedVoronoiCellArea[idx_1];
		double count = 0;

		//double a, b, c;
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
					if ((*fv_it).idx() == idx_1) break;
					index++;
				}
				coe_ij += cot[t][index];
			}
			// Right Triangle
			triangle = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(*voh_it));
			if (!ptr_mesh_->is_boundary(ptr_mesh_->opposite_halfedge_handle(*voh_it)))
			{
				t = triangle.idx();
				index = 0;
				for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(triangle); fv_it.is_valid(); ++fv_it) {
					if ((*fv_it).idx() == idx_1) break;
					index++;
				}
				coe_ij += cot[t][index == 0 ? 2 : index - 1];
			}
			// Set
			count += coe_ij;
			int v_index = ptr_mesh_->to_vertex_handle(*voh_it).idx();
			tripletList.push_back(T(idx_1, v_index, coe_ij / area_2));
		}

		tripletList.push_back(T(idx_1, idx_1, -count / area_2));
	}
		
	Cotangent_Laplacian.setFromTriplets(tripletList.begin(), tripletList.end()); //根据三元组列表生成稀疏矩阵
	Cotangent_Laplacian.makeCompressed();
}

void DifferentialGeometry::Build_BarycentricCell()
{
	current_cell_type_ = BarycentricCell;
	if (is_build_barycentriccell) return;

	BarycentricCellArea.resize(ptr_mesh_->n_vertices());

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		double area = 0;
		MyMesh::VertexVertexIter vbegin = ptr_mesh_->vv_begin(v_it);

		MyMesh::VertexHandle vlast = (*vbegin);
		MyMesh::VertexHandle vbegin_0 = vlast;
		vbegin++;
		Eigen::Vector3d b, c;
		for (; vbegin != ptr_mesh_->vv_end(v_it); vbegin++) {
			auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(*vbegin));
			auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
			b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
			c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
			area += (b.cross(c)).squaredNorm() / 6.0;

			vlast = (*vbegin);
		}
		auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vbegin_0));
		auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
		b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
		c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
		area += (b.cross(c)).squaredNorm() / 6.0;

		BarycentricCellArea[v_it.handle().idx()] = area;
	}

	is_build_barycentriccell = true;
}

void DifferentialGeometry::Build_VoronoiCell()
{
	current_cell_type_ = VoronoiCell;
	if (is_build_voronoicell) return;

	VoronoiCellArea.resize(ptr_mesh_->n_vertices());

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		double area = 0;
		MyMesh::VertexVertexIter vbegin = ptr_mesh_->vv_begin(v_it);

		MyMesh::VertexHandle vlast = (*vbegin);
		MyMesh::VertexHandle vbegin_0 = vlast;
		vbegin++;
		Eigen::Vector3d a, b, c;
		for (; vbegin != ptr_mesh_->vv_end(v_it); vbegin++) {
			auto v_a = (ptr_mesh_->point(*vbegin) - ptr_mesh_->point(vlast));
			auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(*vbegin));
			auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
			a = Eigen::Vector3d(v_a[0], v_a[1], v_a[2]);
			b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
			c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
			double Aa = MY_PI - acos((b.transpose() * c)(0) / (b.norm()*c.norm()));
			double r = a.squaredNorm() / (2.0 * sin(Aa));
			area += (b.cross(c)).squaredNorm() / 4.0 - r * r * sin(2 * Aa) / 4.0;

			vlast = (*vbegin);
		}
		auto v_a = (ptr_mesh_->point(vbegin_0) - ptr_mesh_->point(vlast));
		auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vbegin_0));
		auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
		a = Eigen::Vector3d(v_a[0], v_a[1], v_a[2]);
		b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
		c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
		double Aa = MY_PI - acos((b.transpose() * c)(0) / (b.norm()*c.norm()));
		double r = a.squaredNorm() / (2.0 * sin(Aa));
		area += (b.cross(c)).squaredNorm() / 4.0 - r * r * sin(2 * Aa) / 4.0;

		VoronoiCellArea[(*v_it).idx()] = area;
	}

	is_build_voronoicell = true;
}

void DifferentialGeometry::Build_MixedVoronoiCell()
{
	current_cell_type_ = MixedVoronoiCell;
	if (is_build_mixedvoronoicell) return;
	
	MixedVoronoiCellArea.resize(ptr_mesh_->n_vertices());

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		double area = 0;
		MyMesh::VertexVertexIter vbegin = ptr_mesh_->vv_begin(v_it);

		MyMesh::VertexHandle vlast = (*vbegin);
		MyMesh::VertexHandle vbegin_0 = vlast;
		vbegin++;
		Eigen::Vector3d a, b, c;
		for (; vbegin != ptr_mesh_->vv_end(v_it); vbegin++) {
			auto v_a = (ptr_mesh_->point(*vbegin) - ptr_mesh_->point(vlast));
			auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(*vbegin));
			auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
			a = Eigen::Vector3d(v_a[0], v_a[1], v_a[2]);
			b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
			c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
			if ((a.dot(b) < 0 && b.dot(c) < 0 && c.dot(a) < 0)) {
				double Aa = MY_PI - acos((b.transpose() * c)(0) / (b.norm()*c.norm()));
				double r = a.squaredNorm() / (2.0 * sin(Aa));
				area += (b.cross(c)).squaredNorm() / 4.0 - r * r * sin(2 * Aa) / 4.0;
			}
			else {
				area += (b.cross(c)).squaredNorm() / 4.0;
			}

			vlast = (*vbegin);
		}
		auto v_a = (ptr_mesh_->point(vbegin_0) - ptr_mesh_->point(vlast));
		auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vbegin_0));
		auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(v_it));
		a = Eigen::Vector3d(v_a[0], v_a[1], v_a[2]);
		b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
		c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
		if ((a.dot(b) < 0 && b.dot(c) < 0 && c.dot(a) < 0)) {
			double Aa = MY_PI - acos((b.transpose() * c)(0) / (b.norm()*c.norm()));
			double r = a.squaredNorm() / (2.0 * sin(Aa));
			area += (b.cross(c)).squaredNorm() / 4.0 - r * r * sin(2 * Aa) / 4.0;
		}
		else {
			area += (b.cross(c)).squaredNorm() / 4.0;
		}

		MixedVoronoiCellArea[(*v_it).idx()] = area;
	}
	
	is_build_mixedvoronoicell = true;
}

void DifferentialGeometry::Build_FaceCenter()
{
	FaceCenter.resize(ptr_mesh_->n_faces());
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		MyMesh::Point c = MyMesh::Point(0.0, 0.0, 0.0);
		int count = 0;
		//可能会因为不是逆时针而出问题！！！
		for (auto fv_it = ptr_mesh_->fv_begin(f_it.handle()); fv_it != ptr_mesh_->fv_end(f_it.handle()); ++fv_it) {
			c += ptr_mesh_->point(fv_it);
			count++;
		}
		if (count > 0) c /= (1.0 * count);
		FaceCenter[(*f_it).idx()] = Eigen::Vector3d(c[0],c[1],c[2]);
	}
}

void DifferentialGeometry::Clear_FaceCenter()
{
	FaceCenter.clear();
}

void DifferentialGeometry::Caculate_FaceArea()
{
	average_edgelength = 0.0;
	int count = 0;
	FaceArea.resize(ptr_mesh_->n_faces());
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it)
	{
		double a = 0;
		auto fv_it = ptr_mesh_->fv_begin(f_it.handle());
		auto vbegin = *fv_it;
		fv_it++;
		auto vlast = *fv_it;
		fv_it++;
		Eigen::Vector3d b, c;
		average_edgelength += (ptr_mesh_->point(vlast) - ptr_mesh_->point(vbegin)).length();
		count++;
		//可能会因为不是逆时针而出问题！！！
		for (; fv_it != ptr_mesh_->fv_end(f_it.handle()); ++fv_it) {
			auto v_b = (ptr_mesh_->point(fv_it) - ptr_mesh_->point(vbegin));
			auto v_c = (ptr_mesh_->point(vlast) - ptr_mesh_->point(vbegin));
			b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
			c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
			
			a += (b.cross(c)).norm();
			average_edgelength += (ptr_mesh_->point(fv_it) - ptr_mesh_->point(vlast)).length();
			count++;
			vlast = *fv_it;
		}
		FaceArea[(*f_it).idx()] = a;
	}
	average_edgelength /= count;
}

void DifferentialGeometry::Caculate_Cot()
{
	size_t nT = ptr_mesh_->n_faces();
	cot.resize(nT);
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		int i = (*f_it).idx();
		std::vector<MyMesh::VertexHandle> x_tmp;
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			x_tmp.push_back(*fv_it);
		}

		for (int j = 0; j < 3; j++) {
			auto a = ptr_mesh_->point(x_tmp[j]) - ptr_mesh_->point(x_tmp[j == 0 ? 2 : j - 1]);
			auto b = ptr_mesh_->point(x_tmp[j == 2 ? 0 : j + 1]) - ptr_mesh_->point(x_tmp[j == 0 ? 2 : j - 1]);

			float cos_theta = (a | b) / (a.norm() * b.norm());
			float abs_cot = pow(cos_theta * cos_theta / (1 - cos_theta * cos_theta), 0.5);

			cot[i].push_back(cos_theta > 0 ? abs_cot : -abs_cot);
		}
	}
}

double DifferentialGeometry::Caculate_Volume()
{
	double volume = 0;
	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		Eigen::Matrix3d Vf;
		Vf.setZero();
		int i = 0;
		for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			Vf(i, 0) = (ptr_mesh_->point(*fv_it)[0] - Geometry_Analysis_Center(0));
			Vf(i, 1) = (ptr_mesh_->point(*fv_it)[1] - Geometry_Analysis_Center(1));
			Vf(i, 2) = (ptr_mesh_->point(*fv_it)[2] - Geometry_Analysis_Center(2));
			i++;
		}
		volume += Vf.determinant();
	}
	return volume;
}

void DifferentialGeometry::Caculate_MeanCurvature()
{
	if (is_caculated_meancurvature) return;
	SetCurrentCell(MixedVoronoiCell);
	SetCurrentLaplacian(CotangentLaplacian);

	mean_curvature.resize(ptr_mesh_->n_vertices());
	int row, col;
	row = col = ptr_mesh_->n_vertices();
	Eigen::VectorXd x(row);
	Eigen::VectorXd y(row);
	Eigen::VectorXd z(row);

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		x((*v_it).idx()) = ptr_mesh_->point(v_it)[0];
		y((*v_it).idx()) = ptr_mesh_->point(v_it)[1];
		z((*v_it).idx()) = ptr_mesh_->point(v_it)[2];
	}

	x = Cotangent_Laplacian * x;
	y = Cotangent_Laplacian * y;
	z = Cotangent_Laplacian * z;

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		mean_curvature[(*v_it).idx()] = -0.5 * (x((*v_it).idx()) * ptr_mesh_->normal(*v_it)[0] + y((*v_it).idx()) *  ptr_mesh_->normal(*v_it)[1] + z((*v_it).idx()) *  ptr_mesh_->normal(*v_it)[2]);
	}

	is_caculated_meancurvature = true;
}

void DifferentialGeometry::Caculate_GaussianCurvature()
{
	if (is_caculated_gaussiancurvature) return;
	SetCurrentCell(MixedVoronoiCell);

	gaussian_curvature.resize(ptr_mesh_->n_vertices());

	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); v_it++) {
		double area = MixedVoronoiCellArea[(*v_it).idx()];
		MyMesh::VertexVertexIter vbegin = ptr_mesh_->vv_begin(v_it);

		MyMesh::VertexHandle vlast = (*vbegin);
		MyMesh::VertexHandle vbegin_0 = vlast;
		vbegin++;
		Eigen::Vector3d b, c;
		double theta = 0;
		for (; vbegin != ptr_mesh_->vv_end(v_it); vbegin++) {
			auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(*vbegin));
			auto v_c = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vlast));
			b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
			c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
			
			theta += acos((b.transpose() * c)(0) / (b.norm()*c.norm()));

			vlast = (*vbegin);
		}
		auto v_b = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vbegin_0));
		auto v_c = (ptr_mesh_->point(v_it) - ptr_mesh_->point(vlast));
		b = Eigen::Vector3d(v_b[0], v_b[1], v_b[2]);
		c = Eigen::Vector3d(v_c[0], v_c[1], v_c[2]);
		theta += acos((b.transpose() * c)(0) / (b.norm()*c.norm()));

		gaussian_curvature[(*v_it).idx()] = 1 / area * (2 * MY_PI - theta);
	}

	is_caculated_gaussiancurvature = true;
}

void DifferentialGeometry::Get_MeanCurvature(vector<double> &x)
{
	x.clear();
	for (int i = 0; i < mean_curvature.size(); i++) {
		x.push_back(mean_curvature[i]);
	}
}

void DifferentialGeometry::Get_GaussianCurvature(vector<double>& x)
{
	x.clear();
	for (int i = 0; i < gaussian_curvature.size(); i++) {
		x.push_back(gaussian_curvature[i]);
	}
}

std::tuple<double, double, double> DifferentialGeometry::BarycentricCoordinate(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C, Eigen::Vector3d P)
{
	Eigen::Vector3d n = ((B - A).cross(C - A)).normalized();
	double a = (((B - A).cross(C - A)).transpose()) * n;
	double c1 = (((B - P).cross(C - P)).transpose() * (n / a));
	double c2 = (((C - P).cross(A - P)).transpose() * (n / a));
	double c3 = (((A - P).cross(B - P)).transpose() * (n / a));
	return { c1,c2,c3 };
}
