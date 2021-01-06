#include "Triangulation.h"

Triangulation::Triangulation()
{
}

Triangulation::~Triangulation()
{
}

void Triangulation::triangulate()
{
}

DelaunayTriangulation::DelaunayTriangulation()
{
}

DelaunayTriangulation::DelaunayTriangulation(MyMesh* ptr)
{
	ptr_mesh_ = ptr;
	is_flatten_ = true;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		if (ptr_mesh_->point(v_it)[2] != 0) {
			is_flatten_ = false;
			break;
		}
	}
	if (is_flatten_ && ptr_mesh_->n_edges() > 0) {
		for (MyMesh::EdgeIter e_it = ptr_mesh_->edges_begin(); e_it != ptr_mesh_->edges_end(); ++e_it) {
			ptr_mesh_->delete_edge(e_it, false);
		}
		ptr_mesh_->garbage_collection();
	}
	ptr_mesh_->add_property(is_added);
	/*
	if (!is_flatten_) {
		int fn = ptr_mesh_->n_faces();
		int j = 0;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); j != fn; ++f_it) {
			MyMesh::Point p(0, 0, 0);
			MyMesh::VertexHandle* v = new MyMesh::VertexHandle[3];
			int i = 0;
			for (MyMesh::FaceVertexIter fv_it = ptr_mesh_->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
				v[i++] = *fv_it;
				p+=ptr_mesh_->point(*fv_it);
			}
			p /= 3.0;
			ptr_mesh_->delete_face(f_it, false);
			MyMesh::VertexHandle ph = ptr_mesh_->add_vertex(p);
			ptr_mesh_->add_face(v[0], v[1], ph);
			ptr_mesh_->add_face(v[1], v[2], ph);
			ptr_mesh_->add_face(v[2], v[0], ph);
			j++;
		}
		ptr_mesh_->garbage_collection();
	}
*/
}

DelaunayTriangulation::~DelaunayTriangulation()
{
}

void DelaunayTriangulation::triangulate()
{
	if (!is_flatten_) return;
	init_convexhull();
	init_tri();

	bool s = true;
	while (s) {
		s = false;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); f_it != ptr_mesh_->faces_end(); ++f_it) {
			for (MyMesh::FaceEdgeIter fe_it = ptr_mesh_->fe_begin(*f_it); fe_it != ptr_mesh_->fe_end(*f_it); fe_it++) {
				if (!ptr_mesh_->is_boundary(fe_it)) {
					auto he = ptr_mesh_->halfedge_handle(*fe_it, 0);
					MyMesh::Point v1 = ptr_mesh_->point(ptr_mesh_->from_vertex_handle(he));
					MyMesh::Point v2 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(he));
					MyMesh::Point v0 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(he)));
					MyMesh::Point v3 = ptr_mesh_->point(ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(he))));

					auto e11 = (v1 - v0).normalized();
					auto e12 = (v2 - v0).normalized();
					auto e21 = (v1 - v3).normalized();
					auto e22 = (v2 - v3).normalized();

					float cosa = e11[0] * e12[0] + e11[1] * e12[1];
					float sina = abs(e11[0] * e12[1] - e11[1] * e12[0]);
					float cosb = e21[0] * e22[0] + e21[1] * e22[1];
					float sinb = abs(e22[0] * e21[1] - e22[1] * e21[0]);

					if (cosa * sinb + sina * cosb < 0) {
						if (ptr_mesh_->is_flip_ok(fe_it)) {
							s = true;
							ptr_mesh_->flip(fe_it);
							ptr_mesh_->garbage_collection();
							break;
						}
					}
				}
			}
			if (s) break;
		}
	}
}

bool DelaunayTriangulation::is_inside_tri(MyMesh::VertexHandle p, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::VertexHandle v3, bool strict)
{
	auto p0 = ptr_mesh_->point(p);
	auto p1 = ptr_mesh_->point(v1);
	auto p2 = ptr_mesh_->point(v2);
	auto p3 = ptr_mesh_->point(v3);

	float a1 = (p2 - p1)[0] * (p0 - p1)[1] - (p2 - p1)[1] * (p0 - p1)[0];
	float a2 = (p3 - p2)[0] * (p0 - p2)[1] - (p3 - p2)[1] * (p0 - p2)[0];
	float a3 = (p1 - p3)[0] * (p0 - p3)[1] - (p1 - p3)[1] * (p0 - p3)[0];

	if (strict) {
		if (a1 * a2 > 0 && a2 * a3 > 0) return true;
		else return false;
	}
	else {
		if (a1 * a2 >= 0 && a2 * a3 >= 0 && a3 * a1 >= 0) return true;
		else return false;
	}
}

void DelaunayTriangulation::init_convexhull()
{
	convexhull.clear();
	std::vector<MyMesh::VertexHandle> convexhull_idx;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		bool s = true;
		for (MyMesh::VertexIter v_it1 = ptr_mesh_->vertices_begin(); v_it1 != ptr_mesh_->vertices_end(); ++v_it1) {
			if (v_it1 == v_it) continue;
			for (MyMesh::VertexIter v_it2 = ptr_mesh_->vertices_begin(); v_it2 != ptr_mesh_->vertices_end(); ++v_it2) {
				if (v_it2 == v_it || v_it2 == v_it1) continue;
				for (MyMesh::VertexIter v_it3 = ptr_mesh_->vertices_begin(); v_it3 != ptr_mesh_->vertices_end(); ++v_it3) {
					if (v_it3 == v_it || v_it3 == v_it1 || v_it3 == v_it2) continue;
					if (is_inside_tri(*v_it, *v_it1, *v_it2, *v_it3, true)) {
						s = false;
						break;
					}
				}
				if (!s) break;
			}
			if (!s) break;
		}
		ptr_mesh_->property(is_added, *v_it) = s;
		if(s) convexhull_idx.push_back(*v_it);
	}
	

	while(!convexhull_idx.empty()){
		int idx = convexhull_idx.size() - 1;
		convexhull.push_back(convexhull_idx[idx]);

		if (idx == 0) break;

		for (int j = 0; j < idx; j++) {
			bool s = true;
			std::vector<int> line;
			line.push_back(j);

			auto e1 = ptr_mesh_->point(convexhull_idx[j]) - ptr_mesh_->point(convexhull_idx[idx]);
			for (int k = 0; k < idx; k++) {
				if (k == j) continue;

				auto e2 = ptr_mesh_->point(convexhull_idx[k]) - ptr_mesh_->point(convexhull_idx[idx]);
				float a = e1[0] * e1[1] - e1[1] * e1[0];

				if (abs(a) < 1E-9) {
					line.push_back(k);
				}

				if (a < 0) {
					s = false;
					line.clear();
					break;
				}
			}

			if (s) {
				int idx_next = line[0];
				e1 = ptr_mesh_->point(convexhull_idx[line[0]]) - ptr_mesh_->point(convexhull_idx[idx]);
				float len = e1.length();
				for (int i = 0; i < line.size(); i++) {
					e1 = ptr_mesh_->point(convexhull_idx[line[i]]) - ptr_mesh_->point(convexhull_idx[idx]);
					if (len > e1.length()) {
						idx_next = i;
						len = e1.length();
					}
				}

				convexhull_idx.pop_back();
				if (idx_next != convexhull_idx.size() - 1){
					std::iter_swap(convexhull_idx.begin() + idx_next, convexhull_idx.end() - 1);
				}
				break;
			}
		}
	}
	
	for (int i = 2; i < convexhull.size(); i++) {
		ptr_mesh_->add_face(convexhull[i - 1], convexhull[i], convexhull[0]);
	}
}

void DelaunayTriangulation::init_tri()
{
	std::vector<MyMesh::VertexHandle> add_vertex;
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		if (ptr_mesh_->property(is_added, *v_it)) continue;
		else add_vertex.push_back(*v_it);
	}
	while (!add_vertex.empty()) {
		int fn = ptr_mesh_->n_faces();
		int j = 0;
		for (MyMesh::FaceIter f_it = ptr_mesh_->faces_begin(); j != fn; ++f_it) {
			MyMesh::VertexHandle p = add_vertex.back();
			MyMesh::VertexHandle* v = new MyMesh::VertexHandle[3];
			MyMesh::FaceHalfedgeIter fh_it = ptr_mesh_->fh_iter(*f_it);
			v[0] = ptr_mesh_->from_vertex_handle(fh_it);
			v[1] = ptr_mesh_->to_vertex_handle(fh_it);
			v[2] = ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(fh_it));

			if (is_inside_tri(p, v[0], v[1], v[2], false)) {
				auto e1 = ptr_mesh_->point(v[1]) - ptr_mesh_->point(v[0]);
				auto p1 = ptr_mesh_->point(p) - ptr_mesh_->point(v[0]);
				auto e2 = ptr_mesh_->point(v[2]) - ptr_mesh_->point(v[1]);
				auto p2 = ptr_mesh_->point(p) - ptr_mesh_->point(v[1]);
				auto e3 = ptr_mesh_->point(v[0]) - ptr_mesh_->point(v[2]);
				auto p3 = ptr_mesh_->point(p) - ptr_mesh_->point(v[2]);

				if (abs(e1[0] * p1[1] - e1[1] * p1[0]) > 1E-7) {
					if (abs(e2[0] * p2[1] - e2[1] * p2[0]) > 1E-7) {
						if (abs(e3[0] * p3[1] - e3[1] * p3[0]) > 1E-7) {
							ptr_mesh_->delete_face(f_it, false);
							ptr_mesh_->add_face(v[0], v[1], p);
							ptr_mesh_->add_face(v[1], v[2], p);
							ptr_mesh_->add_face(v[2], v[0], p);
						}
						else {
							auto fh_it_ = ptr_mesh_->next_halfedge_handle(ptr_mesh_->next_halfedge_handle(fh_it));
							auto ff_it = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it_));
							auto v_ = ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(fh_it_)));

							ptr_mesh_->delete_face(f_it, false);
							ptr_mesh_->delete_face(ff_it, false);
							ptr_mesh_->add_face(v_, v[0], p);
							ptr_mesh_->add_face(v[0], v[1], p);
							ptr_mesh_->add_face(v[1], v[2], p);
							ptr_mesh_->add_face(v[2], v_, p);
						}						
					}
					else {
						auto fh_it_ = ptr_mesh_->next_halfedge_handle(fh_it);
						auto ff_it = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it_));
						auto v_ = ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(fh_it_)));

						ptr_mesh_->delete_face(f_it, false);
						ptr_mesh_->delete_face(ff_it, false);
						ptr_mesh_->add_face(v_, v[2], p);
						ptr_mesh_->add_face(v[2], v[0], p);
						ptr_mesh_->add_face(v[0], v[1], p);
						ptr_mesh_->add_face(v[1], v_, p);
					}
				}
				else {
					auto ff_it = ptr_mesh_->face_handle(ptr_mesh_->opposite_halfedge_handle(fh_it));
					auto v_ = ptr_mesh_->to_vertex_handle(ptr_mesh_->next_halfedge_handle(ptr_mesh_->opposite_halfedge_handle(fh_it)));
					
					ptr_mesh_->delete_face(f_it, false);
					ptr_mesh_->delete_face(ff_it, false);
					ptr_mesh_->add_face(v_, v[1], p);
					ptr_mesh_->add_face(v[1], v[2], p);
					ptr_mesh_->add_face(v[2], v[0], p);
					ptr_mesh_->add_face(v[0], v_, p);
				}
				
				add_vertex.pop_back();
				break;
			}
			j++;
		}
		if (j == fn) {
			auto v = add_vertex.back();
			add_vertex.pop_back();
			add_vertex.insert(add_vertex.begin(), v);
		}
		ptr_mesh_->garbage_collection();
	}
}
