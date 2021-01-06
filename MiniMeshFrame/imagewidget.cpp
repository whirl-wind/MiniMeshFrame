#include "imagewidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include "Parameterization.h"

using std::cout;
using std::endl;

ImageWidget::ImageWidget(ImageWindow *mainwindow)
{
	ptr_imagewindow_ = mainwindow;
	ptr_mesh_ = ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->ptr_mesh_;
	if(Open(ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->texture_filename_)) image_status_ = true;
	Pt_Uniform();
}


ImageWidget::~ImageWidget(void)
{
}

void ImageWidget::paintEvent(QPaintEvent *paintevent)
{
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);
	int ax = (width() - 250) / 2;
	int ay = (height() - 250) / 2;
	// Draw image
	if (image_status_) {
		int w = image_mat_.cols, h = image_mat_.rows;
		if (w > width()) {
			h = width()*h / w;
			w = width();
		}
		if (h > height()) {
			w = height()*w / h;
			h = height();
		}
		QImage image_show = QImage((unsigned char *)(image_mat_.data), image_mat_.cols, image_mat_.rows, image_mat_.step, QImage::Format_RGB888);
		QImage temp = image_show.scaled(w, h, Qt::KeepAspectRatio);
		QRect rect = QRect((width() - temp.width()) / 2, (height() - temp.height()) / 2, temp.width(), temp.height());
		painter.drawImage(rect, image_show);
		ax = (width() - temp.width()) / 2;
		ay = (height() - temp.height()) / 2;
	}

	QColor green(0, 255, 0);
	QColor blue(0, 0, 255);
	QColor greenblue(0, 191, 191, 191);
	QColor red(191, 0, 0, 191);
	QPen pen1(green);
	QPen pen2(blue);
	QPen pen3(greenblue);
	QPen pen4(red);
	pen4.setWidth(2);
	pen3.setWidth(2);
	pen2.setWidth(5);
	pen1.setWidth(5);

	painter.setPen(pen3);
	painter.drawLine(start_point_, end_point_);
	painter.setPen(pen2);
	painter.drawPoint(start_point_);
	painter.setPen(pen1);
	painter.drawPoint(end_point_);

	pt->Draw(painter,ax,ay);

	update();
	painter.end();
}

void ImageWidget::Open()
{
	// Open file
	filename_ = QFileDialog::getOpenFileName(this, tr("Read Image"), "D:\\Zz\\CG\\Materials", tr("Images(*.bmp *.png *.jpg)"));

	// Load file
	if (!filename_.isEmpty())
	{
		image_mat_ = cv::imread(filename_.toLatin1().data());
		cvtColor(image_mat_, image_mat_, cv::COLOR_BGR2RGB);
		image_mat_backup_ = image_mat_.clone();
		image_status_ = true;
		pt = new Pt_uniform(ptr_mesh_, image_mat_.cols, image_mat_.rows);
		pt->parameterization();
	}
	else
		return;
	update();
}

void ImageWidget::Open(Mat img)
{
	image_mat_ = img;
	image_mat_backup_ = image_mat_.clone();
	image_status_ = true;
	pt = new Pt_uniform(ptr_mesh_, image_mat_.cols, image_mat_.rows);
	pt->parameterization();

	update();
}

bool ImageWidget::Open(QString img_filename_)
{
	// Open file
	filename_ = img_filename_;

	// Load file
	if (!filename_.isEmpty())
	{
		image_mat_ = cv::imread(filename_.toLatin1().data());
		cvtColor(image_mat_, image_mat_, cv::COLOR_BGR2RGB);
		image_mat_backup_ = image_mat_.clone();
		image_status_ = true;
		pt = new Pt_uniform(ptr_mesh_, image_mat_.cols, image_mat_.rows);
		pt->parameterization();
	}
	else
		return false;
	update();
	return true;
}

void ImageWidget::Save()
{
	if (filename_.isNull())
		SaveAs();
	else {
		Mat image_save;
		cvtColor(image_mat_, image_save, cv::COLOR_BGR2RGB);
		imwrite(filename_.toLatin1().data(), image_save);
	}
}

void ImageWidget::SaveAs()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), ".", tr("Images(*.bmp *.png *.jpg)"));
	if (filename.isNull())
	{
		return;
	}
	QImage image_show = QImage((unsigned char *)(image_mat_.data), image_mat_.cols, image_mat_.rows, image_mat_.step, QImage::Format_RGB888);
	QPixmap pixmap(image_show.width()+1, image_show.height()+1);
	this->render(&pixmap, QPoint(), QRegion((width() - image_show.width()) / 2, (height() - image_show.height()) / 2, image_show.width()+1, image_show.height()+1));//保存当前窗口
	pixmap.save(filename);
}

void ImageWidget::mousePressEvent(QMouseEvent * event)
{
	if (Qt::LeftButton == event->button()) {
		start_point_ = end_point_ = event->pos();
	}
}

void ImageWidget::mouseMoveEvent(QMouseEvent * event)
{
	end_point_ = event->pos();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent * event)
{

}

void ImageWidget::Invert()
{
	MatIterator_<Vec3b> iter, iterend;
	for (iter = image_mat_.begin<Vec3b>(), iterend = image_mat_.end<Vec3b>(); iter != iterend; ++iter)
	{
		(*iter)[0] = 255 - (*iter)[0];
		(*iter)[1] = 255 - (*iter)[1];
		(*iter)[2] = 255 - (*iter)[2];
	}
	update();
}

void ImageWidget::Restore()
{
	image_mat_ = image_mat_backup_.clone();
	update();
}

void ImageWidget::Pt_Uniform()
{
	pt = new Pt_uniform(ptr_mesh_, 250, 250);
	if(image_status_) pt = new Pt_uniform(ptr_mesh_, image_mat_.cols, image_mat_.rows);
	if (!pt->is_good_) return;
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}

void ImageWidget::Pt_Weightedleastsquares()
{
	pt = new Pt_weightedleastsquares(ptr_mesh_, 250, 250);
	if (image_status_) pt = new Pt_weightedleastsquares(ptr_mesh_, image_mat_.cols, image_mat_.rows);
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}

void ImageWidget::Pt_Shapepreserving()
{
	pt = new Pt_shapepreserving(ptr_mesh_, 250, 250);
	if (image_status_) pt = new Pt_shapepreserving(ptr_mesh_, image_mat_.cols, image_mat_.rows);
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}

void ImageWidget::Pt_Asap()
{
	pt = new Pt_ASAP(ptr_mesh_, 250, 250, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	if (image_status_) pt = new Pt_ASAP(ptr_mesh_, image_mat_.cols, image_mat_.rows, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}

void ImageWidget::Pt_Arap()
{
	pt = new Pt_ARAP(ptr_mesh_, 250, 250, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	if (image_status_) pt = new Pt_ARAP(ptr_mesh_, image_mat_.cols, image_mat_.rows, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}

void ImageWidget::Pt_Slim()
{
	pt = new Pt_SLIM(ptr_mesh_, 250, 250, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	if (image_status_) pt = new Pt_SLIM(ptr_mesh_, image_mat_.cols, image_mat_.rows, this->ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->DDG);
	pt->parameterization();
	if (image_status_) pt->SetTexcoord();
	update();
	ptr_imagewindow_->ptr_mainwindow_->ptr_renderingwidget_->updateGL();
}
