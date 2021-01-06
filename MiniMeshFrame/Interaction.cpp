#include "Interaction.h"
#include <gl/glut.h>
#include <gl/GLU.h>

Interaction::Interaction()
{
	point_array.clear();
	center = MyMesh::Point(0,0,0);
	fix = false;
	selected_direction = 3;
}

Interaction::Interaction(MyMesh* ptr_object_, bool is_fix)
{
	ptr_mesh_ = ptr_object_;
	float x, y, z;
	x = y = z = 0.0;
	point_array.clear();
	for (MyMesh::VertexIter v_it = ptr_mesh_->vertices_begin(); v_it != ptr_mesh_->vertices_end(); ++v_it) {
		point_array.push_back(*v_it);
		x += ptr_mesh_->point(*v_it)[0];
		y += ptr_mesh_->point(*v_it)[1];
		z += ptr_mesh_->point(*v_it)[2];
	}
	float n = ptr_mesh_->n_vertices();
	center = MyMesh::Point(x / n, y / n, z / n);

	fix = is_fix;
	draw_point_ = false;
	selected_direction = 3;
}

Interaction::Interaction(MyMesh* ptr_object_, vector<MyMesh::VertexHandle> p_a, bool is_fix)
{
	ptr_mesh_ = ptr_object_;
	float x, y, z;
	x = y = z = 0.0;
	point_array.clear();
	//point_array.assign(p_a.begin(), p_a.end());
	float n = 0;
	for (vector<MyMesh::VertexHandle>::iterator v_it = p_a.begin(); v_it != p_a.end(); ++v_it) {
		point_array.push_back(*v_it);
		x += ptr_mesh_->point(*v_it)[0];
		y += ptr_mesh_->point(*v_it)[1];
		z += ptr_mesh_->point(*v_it)[2];
		n += 1.0;
	}
	center = MyMesh::Point(x / n, y / n, z / n);

	fix = is_fix;
	draw_point_ = true;
	selected_direction = 3;
}

Interaction::~Interaction()
{
}

void Interaction::update_selection()
{
	selected_direction = 3;
	auto ind = handle_point - center;
	for (int i = 0; i < 3; i++) {
		if (ind[i]>0 && (abs(ind[(i + 1) % 3]) + abs(ind[(i + 2) % 3]) <= abs(ind[(selected_direction + 1) % 3]) + abs(ind[(selected_direction + 2) % 3]))) selected_direction = i;
	}
}

void Interaction::update_center()
{
	float x, y, z;
	x = y = z = 0.0;
	//point_array.assign(p_a.begin(), p_a.end());
	float n = 0;
	for (vector<MyMesh::VertexHandle>::iterator v_it = point_array.begin(); v_it != point_array.end(); ++v_it) {
		x += ptr_mesh_->point(*v_it)[0];
		y += ptr_mesh_->point(*v_it)[1];
		z += ptr_mesh_->point(*v_it)[2];
		n += 1.0;
	}
	center = MyMesh::Point(x / n, y / n, z / n);
}

void Interaction::Draw()
{
	if (!fix) {
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex3f(center[0], center[1], center[2]);
		glVertex3f(center[0] + 0.5, center[1], center[2]);
		glEnd();
		glPushMatrix();
		glTranslatef(center[0] + 0.5, center[1], center[2]);
		glRotatef(90, 0.0, 1.0, 0.0);
		//glutSolidCone(0.02, 0.06, 20, 10);
		glPopMatrix();

		//y axis
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3f(center[0], center[1], center[2]);
		glVertex3f(center[0], center[1] + 0.5, center[2]);
		glEnd();
		glPushMatrix();
		glTranslatef(center[0], center[1] + 0.5, center[2]);
		glRotatef(90, -1.0, 0.0, 0.0);
		//glutSolidCone(0.02, 0.06, 20, 10);
		glPopMatrix();

		//z axis
		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);
		glVertex3f(center[0], center[1], center[2]);
		glVertex3f(center[0], center[1], center[2] + 0.5);
		glEnd();
		glPushMatrix();
		glTranslatef(center[0], center[1], center[2] + 0.5);
		//glutSolidCone(0.02, 0.06, 20, 10);
		glPopMatrix();
	}

	if (draw_point_) {
		if (fix) glColor3f(1.0, 0.5, 0.0);
		else glColor3f(0.5, 1.0, 0.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < point_array.size(); i++) {
			glNormal3fv(ptr_mesh_->normal(point_array[i]).data());
			glVertex3fv(ptr_mesh_->point(point_array[i]).data());
		}
		glEnd();
	}

	glColor3f(1.0, 1.0, 1.0);
}

void Interaction::Draw_Select()
{
	/*
	glColor3f(1.0, 0.5, 0);
	glPointSize(3);
	glBegin(GL_POINTS);
	glVertex3fv(handle_point.data());
	glEnd();
	glPointSize(1);
	*/
	if (!fix) {
		glLineWidth(2);
		switch (selected_direction) {
		case 0:
			glColor3f(1.0, 0.5, 0.5);
			glBegin(GL_LINES);
			glVertex3f(center[0], center[1], center[2]);
			glVertex3f(center[0] + 0.5, center[1], center[2]);
			glEnd();
			glPushMatrix();
			glTranslatef(center[0] + 0.5, center[1], center[2]);
			glRotatef(90, 0.0, 1.0, 0.0);
			//glutSolidCone(0.02, 0.06, 20, 10);
			glPopMatrix();
			break;
		case 1:
			glColor3f(0.5, 1.0, 0.5);
			glBegin(GL_LINES);
			glVertex3f(center[0], center[1], center[2]);
			glVertex3f(center[0], center[1] + 0.5, center[2]);
			glEnd();
			glPushMatrix();
			glTranslatef(center[0], center[1] + 0.5, center[2]);
			glRotatef(90, -1.0, 0.0, 0.0);
			//glutSolidCone(0.02, 0.06, 20, 10);
			glPopMatrix();
			break;
		case 2:
			glColor3f(0.5, 0.5, 1.0);
			glBegin(GL_LINES);
			glVertex3f(center[0], center[1], center[2]);
			glVertex3f(center[0], center[1], center[2] + 0.5);
			glEnd();
			glPushMatrix();
			glTranslatef(center[0], center[1], center[2] + 0.5);
			//glutSolidCone(0.02, 0.06, 20, 10);
			glPopMatrix();
			break;
		default:
			break;
		}
		glLineWidth(1);
	}
	
	if (draw_point_) {
		if (fix) glColor3f(1.0, 0.25, 0.25);
		else glColor3f(0.25, 1.0, 0.25);
		glPointSize(3);
		glBegin(GL_POINTS);
		for (int i = 0; i < point_array.size(); i++) {
			glNormal3fv(ptr_mesh_->normal(point_array[i]).data());
			glVertex3fv(ptr_mesh_->point(point_array[i]).data());
		}
		glEnd();
		glPointSize(1);
	}

	glColor3f(1.0, 1.0, 1.0);
}

void Interaction::Move(MyMesh::Point d)
{
	if (fix) return;
	if (selected_direction == 3) return;
	MyMesh::Point dn;
	switch (selected_direction) {
	case 0:
		dn = MyMesh::Point(d[0], 0, 0);
		break;
	case 1:
		dn = MyMesh::Point(0, d[1], 0);
		break;
	case 2:
		dn = MyMesh::Point(0, 0, d[2]);
		break;
	}
	for (int i = 0; i < point_array.size(); i++) {
		MyMesh::Point new_point = ptr_mesh_->point(point_array[i]) + dn;
		ptr_mesh_->set_point(point_array[i], new_point);
	}
	center += dn;
}
