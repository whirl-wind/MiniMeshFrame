#include "Remeshing.h"

Remeshing::Remeshing()
{
}

Remeshing::~Remeshing()
{
}

void Remeshing::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Remeshing::remesh()
{
}

IncrementalRemeshing::IncrementalRemeshing(MyMesh* ptr, DifferentialGeometry* ddg)
{
	ptr_mesh_ = ptr;
	LoadDDG(ddg);
}

IncrementalRemeshing::~IncrementalRemeshing()
{
}

void IncrementalRemeshing::remesh(double target_edge_length)
{
	double low = 4.0 / 5.0 * target_edge_length;
	double high = 4.0 / 3.0 * target_edge_length;
	for (size_t i = 0; i < 10; i++) {
		split_long_edges(high);
		collapse_short_edges(low, high);
		equalize_valences();
		tangential_relaxation();
		project_to_surface();
		std::cout << i << std::endl;
	}
	std::cout << "Done\n" << std::endl;
}

void IncrementalRemeshing::split_long_edges(double high)
{
	for (MyMesh::EdgeIter e_it = ptr_mesh_->edges_begin(); e_it != ptr_mesh_->edges_end(); ++e_it) {
		auto he_it = ptr_mesh_->halfedge_handle(e_it, 0);
		MyMesh::Point v1 = ptr_mesh_->point(ptr_mesh_->from_vertex_handle(he_it));
		MyMesh::Point v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(he_it));
		if ((v1 - v2).norm() > high) {
			MyMesh::Point v = (v1 + v2) / 2.0;
			MyMesh::VertexHandle v_h = ptr_mesh_->add_vertex(v);
			ptr_mesh_->split_edge(e_it, v_h);
		}
	}
	ptr_mesh_->garbage_collection();
	std::cout << "split Done\n" << std::endl;
}

void IncrementalRemeshing::collapse_short_edges(double low, double high)
{
	for (MyMesh::EdgeIter e_it = ptr_mesh_->edges_begin(); e_it != ptr_mesh_->edges_end(); ++e_it) {
		//glNormal3fv(ptr_mesh_->normal(ptr_mesh_->from_vertex_handle(*he_it)).data());
		auto he_it = ptr_mesh_->halfedge_handle(e_it, 0);
		MyMesh::VertexHandle a = ptr_mesh_->from_vertex_handle(he_it);
		MyMesh::Point v1 = ptr_mesh_->point(a);
		MyMesh::Point v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(he_it));
		if ((v1 - v2).norm() < low) {
			bool collapse_ok = ptr_mesh_->is_collapse_ok(he_it);
			for (MyMesh::VertexVertexIter vv_it = ptr_mesh_->vv_begin(a); vv_it != ptr_mesh_->vv_end(a); vv_it++) {
				if ((ptr_mesh_->point(vv_it) - v2).norm() > high) collapse_ok = false;
			}
			if (collapse_ok) ptr_mesh_->collapse(he_it);
		}
	}
	ptr_mesh_->garbage_collection();
	std::cout << "collapse Done\n" << std::endl;
}

void IncrementalRemeshing::equalize_valences()
{
	for (MyMesh::EdgeIter e_it = ptr_mesh_->edges_begin(); e_it != ptr_mesh_->edges_end(); ++e_it) {
		if (!ptr_mesh_->is_boundary(*e_it) && (ptr_mesh_->is_flip_ok(*e_it))) {
			auto he_it = ptr_mesh_->halfedge_handle(e_it, 0);
			MyMesh::VertexHandle v1 = ptr_mesh_->from_vertex_handle(he_it);
			MyMesh::VertexHandle v2 = ptr_mesh_->to_vertex_handle(he_it);
			auto next_he = ptr_mesh_->next_halfedge_handle(he_it);
			MyMesh::VertexHandle v3 = ptr_mesh_->to_vertex_handle(next_he);
			next_he = ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(he_it));
			MyMesh::VertexHandle v4 = ptr_mesh_->to_vertex_handle(next_he);

			int vv1, vv2, vv3, vv4;
			if (ptr_mesh_->is_boundary(v1)) vv1 = 4;
			else vv1 = 6;
			if (ptr_mesh_->is_boundary(v2)) vv2 = 4;
			else vv2 = 6;
			if (ptr_mesh_->is_boundary(v3)) vv3 = 4;
			else vv3 = 6;
			if (ptr_mesh_->is_boundary(v4)) vv4 = 4;
			else vv4 = 6;
			
			int deviation_pre =
				  abs(int(ptr_mesh_->valence(v1)) - vv1)
				+ abs(int(ptr_mesh_->valence(v2)) - vv2)
				+ abs(int(ptr_mesh_->valence(v3)) - vv3)
				+ abs(int(ptr_mesh_->valence(v4)) - vv4);

			// Flip edge
			ptr_mesh_->flip(*e_it);
			
			int deviation_post =
				  abs(int(ptr_mesh_->valence(v1)) - vv1)
				+ abs(int(ptr_mesh_->valence(v2)) - vv2)
				+ abs(int(ptr_mesh_->valence(v3)) - vv3)
				+ abs(int(ptr_mesh_->valence(v4)) - vv4);

			//std::cout << deviation_pre << " " << deviation_post << std::endl;
			if (deviation_pre <= deviation_post) {
				ptr_mesh_->flip(*e_it);
			}
		}
		//std::cout <<"valence"<< std::endl;
	}
	ptr_mesh_->garbage_collection();
	std::cout << "valence Done\n" << std::endl;
}

void IncrementalRemeshing::tangential_relaxation()
{
	ptr_mesh_->update_normals();
	std::vector<MyMesh::Point> q;
	q.resize(ptr_mesh_->n_vertices());
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		MyMesh::Point vq(0, 0, 0);
		int count = 0;
		for (MyMesh::VertexVertexIter vv_it = ptr_mesh_->vv_begin(v_it); vv_it != ptr_mesh_->vv_end(v_it); vv_it++) {
			vq += ptr_mesh_->point(vv_it);
			count++;
		}
		vq /= float(count);
		q[(*v_it).idx()] = vq;
	}
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		Eigen::Vector3d q_(q[(*v_it).idx()][0], q[(*v_it).idx()][1], q[(*v_it).idx()][2]);
		Eigen::Vector3d p(ptr_mesh_->point(v_it)[0], ptr_mesh_->point(v_it)[1], ptr_mesh_->point(v_it)[2]);
		Eigen::Vector3d n(ptr_mesh_->normal(v_it)[0], ptr_mesh_->normal(v_it)[1], ptr_mesh_->normal(v_it)[2]);
		Eigen::Vector3d p_ = q_ + (n * n.transpose()) * (p - q_);
		ptr_mesh_->set_point(v_it, MyMesh::Point(p_(0), p_(1), p_(2)));
	}
	std::cout << "tangential Done\n" << std::endl;
}

void IncrementalRemeshing::project_to_surface()
{
}
