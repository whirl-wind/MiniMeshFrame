#pragma once
#ifndef IMAGEWIDGET_H
#define IMAGEWIDGET_H

#include <QWidget>
#include "imagewindow.h"
#include "renderingwidget.h"
#include "Parameterization.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

QT_BEGIN_NAMESPACE
class QImage;
class QPainter;
QT_END_NAMESPACE

using namespace cv;
using std::vector;

enum WarpType
{
	None = -1,
	DefaultW = 0,
	Idw = 1,
	Rbf = 2,
};

class ImageWidget :
	public QWidget
{
	Q_OBJECT

public:
	ImageWidget(ImageWindow *mainwindow = 0);
	~ImageWidget(void);

protected:
	void paintEvent(QPaintEvent *paintevent);

public slots:
	// File IO
	void Open();
	void Open(Mat img);
	bool Open(QString img_filename_);
	// Open an image file, support ".bmp, .png, .jpg" format
	void Save();												// Save image to current file
	void SaveAs();												// Save image to another file

	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);

	// Image processing
	void Invert();												// Invert pixel value in image
	void Restore();												// Restore image to origin

	void Pt_Uniform();
	void Pt_Weightedleastsquares();
	void Pt_Shapepreserving();
	void Pt_Asap();
	void Pt_Arap();
	void Pt_Slim();

private:
	Parameterization *pt;

	Mat			image_mat_;
	Mat			image_mat_backup_;

	QPoint		start_point_;
	QPoint		end_point_;
	//Line* current_line_;

	QString		filename_;
	bool		image_status_;

public:
	ImageWindow	*ptr_imagewindow_;
	MyMesh	*ptr_mesh_;
};

#endif