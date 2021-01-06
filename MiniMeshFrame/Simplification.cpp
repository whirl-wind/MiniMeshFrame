#include "Simplification.h"

Simplification::Simplification()
{
	ptr_mesh_ = new MyMesh();
	DDG = new DifferentialGeometry();
}

Simplification::~Simplification()
{
}

void Simplification::LoadDDG(DifferentialGeometry* ddg)
{
	DDG = ddg;
}

void Simplification::simplify()
{
}

void Simplification::simplify_one_step()
{
}

Sp_QEM::Sp_QEM(MyMesh* ptr_, int n, float t)
{
	ptr_mesh_ = ptr_;
	Nv = n;
	Cboundary = t;
	ptr_mesh_->add_property(QEM_f);
	ptr_mesh_->add_property(QEM_v);
	ptr_mesh_->add_property(QEM_e);
	ptr_mesh_->add_property(ev);

	Init();
}

Sp_QEM::~Sp_QEM()
{
	ptr_mesh_->remove_property(QEM_f);
	ptr_mesh_->remove_property(QEM_v);
	ptr_mesh_->remove_property(QEM_e);
	ptr_mesh_->remove_property(ev);
}

void Sp_QEM::Init()
{
	while (!q.empty()) q.pop();
	N_collapse = 0;
	if (ptr_mesh_ == NULL) return;
	if (ptr_mesh_->n_faces() == 0) return;

	for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
		Eigen::Vector3f n = Eigen::Vector3f(ptr_mesh_->normal(*f_it).data());
		auto vv = ptr_mesh_->fv_begin(f_it).handle();
		Eigen::Vector3f v = Eigen::Vector3f(ptr_mesh_->point(vv).data());
		float d = -1.0 * (n.transpose() * v)(0);
		Eigen::Vector4f N(n(0),n(1),n(2),d);
		ptr_mesh_->property(QEM_f, *f_it) = N * N.transpose();
		//std::cout << N.transpose() << std::endl;
	}
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		Eigen::Matrix4f Q;
		Q.setZero();
		for (MyMesh::VertexFaceIter vf_it = ptr_mesh_->vf_begin(v_it); vf_it != ptr_mesh_->vf_end(v_it); ++vf_it) {
			Q += ptr_mesh_->property(QEM_f, *vf_it);
		}
		ptr_mesh_->property(QEM_v, *v_it) = Q;
	}
	for (MyMesh::EdgeIter e_it = ptr_mesh_->edges_begin(); e_it != ptr_mesh_->edges_end(); ++e_it) {
		auto he_it = ptr_mesh_->halfedge_handle(*e_it, 0);
		auto v1 = ptr_mesh_->from_vertex_handle(he_it);
		auto v2 = ptr_mesh_->to_vertex_handle(he_it);
		ptr_mesh_->property(QEM_e, *e_it) = ptr_mesh_->property(QEM_v, v1) + ptr_mesh_->property(QEM_v, v2);
		
		Eigen::Matrix4f Q = ptr_mesh_->property(QEM_e, *e_it);
		Q(3, 0) = 0;
		Q(3, 1) = 0;
		Q(3, 2) = 0;
		Q(3, 3) = 1;
		Eigen::Vector4f b = Eigen::Vector4f(0, 0, 0, 1);
		Eigen::Vector4f v = Q.lu().solve(b);
		float cost = v.transpose() * ptr_mesh_->property(QEM_e, *e_it) * v;
		//std::cout << cost <<std::endl;
		edge_v* node = new edge_v(cost, *e_it, v);
		q.push(node);
		ptr_mesh_->property(ev, *e_it) = node;
	}
	//while (!q.empty()) printf("(%d,%f)--(%d,%f)\n", q.top()->edge.idx(), q.top()->cost,ptr_mesh_->edge_handle(q.top()->edge.idx()).idx(), ptr_mesh_->property(ev, ptr_mesh_->edge_handle(q.top()->edge.idx()))->cost), q.pop();
}

bool Sp_QEM::Qualified()
{
	float cost = update_min_cost();
	return !q.empty() && cost < Cboundary && Nv + N_collapse < ptr_mesh_->n_vertices();
}

float Sp_QEM::update_min_cost()
{
	while (q.top()->cost == INT_MAX) {
		auto p = q.top();
		q.pop();
		delete p;
		if (q.empty()) return INT_MAX;
	}
	return q.top()->cost;
}

void Sp_QEM::simplify()
{
	while (Qualified())
	{
		simplify_one_step();
	}
	ptr_mesh_->garbage_collection();
	ptr_mesh_->update_normals();
	Init();
}

void Sp_QEM::simplify_one_step()
{
	MyMesh::EdgeHandle e_it = q.top()->edge;
	auto he_it = ptr_mesh_->halfedge_handle(e_it, 0);
		
	vector<edge_v*> nodes;
	while (!ptr_mesh_->is_collapse_ok(he_it)) {
		auto node = q.top();
		nodes.push_back(node);
		q.pop();
		if (!q.empty()) {
			e_it = q.top()->edge;
			he_it = ptr_mesh_->halfedge_handle(e_it, 0);
		}
		else return;
	}
	auto v1 = ptr_mesh_->from_vertex_handle(he_it);
	auto v2 = ptr_mesh_->to_vertex_handle(he_it);
	for (MyMesh::VertexEdgeIter ve_it = ptr_mesh_->ve_begin(v1); ve_it != ptr_mesh_->ve_end(v1); ++ve_it) {
		ptr_mesh_->property(ev, *ve_it)->cost = INT_MAX;
	}
	for (MyMesh::VertexEdgeIter ve_it = ptr_mesh_->ve_begin(v2); ve_it != ptr_mesh_->ve_end(v2); ++ve_it) {
		ptr_mesh_->property(ev, *ve_it)->cost = INT_MAX;
	}
	
	ptr_mesh_->set_point(v2, MyMesh::Point(q.top()->v(0), q.top()->v(1), q.top()->v(2)));
	ptr_mesh_->property(QEM_v, v2) = ptr_mesh_->property(QEM_e, e_it);
	ptr_mesh_->collapse(he_it);
	//ptr_mesh_->delete_vertex(v1);
	N_collapse++;
	q.pop();

	for (MyMesh::VertexFaceIter vf_it = ptr_mesh_->vf_begin(v2); vf_it != ptr_mesh_->vf_end(v2); ++vf_it) {
		MyMesh::HalfedgeHandle h;
		for (MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_begin(vf_it); fh_it != ptr_mesh_->fh_end(vf_it); ++fh_it) {
			h = *fh_it;
			if (ptr_mesh_->to_vertex_handle(h) == v2) break;
		}
		Eigen::Vector3f e1 = Eigen::Vector3f((ptr_mesh_->point(ptr_mesh_->from_vertex_handle(h)) - ptr_mesh_->point(v2)).data());
		Eigen::Vector3f e2 = Eigen::Vector3f((ptr_mesh_->point(ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(h))) - ptr_mesh_->point(v2)).data());
		Eigen::Vector3f n = e1.cross(e2).normalized();
		Eigen::Vector3f v = Eigen::Vector3f(ptr_mesh_->point(v2).data());
		float d = -1.0 * (n.transpose() * v)(0);
		Eigen::Vector4f N(n(0), n(1), n(2), d);
		ptr_mesh_->property(QEM_f, *vf_it) = N * N.transpose();
	}
	Eigen::Matrix4f Q2;
	Q2.setZero();
	for (MyMesh::VertexFaceIter vf_it = ptr_mesh_->vf_begin(v2); vf_it != ptr_mesh_->vf_end(v2); ++vf_it) {
		Q2 += ptr_mesh_->property(QEM_f, *vf_it);
	}
	ptr_mesh_->property(QEM_v, v2) = Q2;
	for (MyMesh::VertexVertexIter vv_it = ptr_mesh_->vv_begin(v2); vv_it != ptr_mesh_->vv_end(v2); ++vv_it) {
		Eigen::Matrix4f Q;
		Q.setZero();
		for (MyMesh::VertexFaceIter vf_it = ptr_mesh_->vf_begin(vv_it); vf_it != ptr_mesh_->vf_end(vv_it); ++vf_it) {
			Q += ptr_mesh_->property(QEM_f, *vf_it);
		}
		ptr_mesh_->property(QEM_v, *vv_it) = Q;
	}
	for (MyMesh::VertexEdgeIter ve_it = ptr_mesh_->ve_begin(v2); ve_it != ptr_mesh_->ve_end(v2); ++ve_it) {
		auto h_it = ptr_mesh_->halfedge_handle(ve_it, 0);
		auto v_1 = ptr_mesh_->from_vertex_handle(h_it);
		auto v_2 = ptr_mesh_->to_vertex_handle(h_it);
		ptr_mesh_->property(QEM_e, ve_it) = ptr_mesh_->property(QEM_v, v_1) + ptr_mesh_->property(QEM_v, v_2);

		Eigen::Matrix4f Q = ptr_mesh_->property(QEM_e, ve_it);
		Q(3, 0) = 0;
		Q(3, 1) = 0;
		Q(3, 2) = 0;
		Q(3, 3) = 1;
		Eigen::Vector4f b = Eigen::Vector4f(0, 0, 0, 1);
		Eigen::Vector4f v = Q.lu().solve(b);
		float cost = v.transpose() * ptr_mesh_->property(QEM_e, ve_it) * v;

		edge_v* node = new edge_v(cost, *ve_it, v);
		q.push(node);
	}
	
	while (!nodes.empty()) {
		auto node = nodes.back();
		q.push(node);
		nodes.pop_back();
	}
}
