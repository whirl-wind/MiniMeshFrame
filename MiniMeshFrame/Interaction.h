#pragma once
#include <QObject>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
using std::vector;

class Interaction : public QObject
{
public:
	Interaction();
	Interaction(MyMesh* ptr_object_, bool is_fix = false);
	Interaction(MyMesh* ptr_object_, vector<MyMesh::VertexHandle> p_a, bool is_fix = false);
	~Interaction();
	
	void update_selection();
	void update_center();
	void Draw();
	void Draw_Select();
	void Move(MyMesh::Point d);

	MyMesh::Point center;
	MyMesh* ptr_mesh_;
	MyMesh::Point handle_point;
	MyMesh::Point near_point, far_point;

	bool fix;
	vector<MyMesh::VertexHandle> point_array;

private:	
	bool draw_point_;
	int selected_direction;
};

