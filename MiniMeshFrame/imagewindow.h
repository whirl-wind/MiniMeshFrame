#pragma once
#ifndef IMAGEWINDOW_H
#define IMAGEWINDOW_H

#include "mainwindow.h"
#include <QtWidgets/QMainWindow>

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class ViewWidget;
class QImage;
class QPainter;
class QRect;
class ImageWidget;
QT_END_NAMESPACE

class ImageWindow : public QMainWindow
{
	Q_OBJECT

public:
	ImageWindow(MainWindow *mainwindow = 0);
	~ImageWindow();

	MainWindow	*ptr_mainwindow_;
	ImageWidget	*imagewidget_;

protected:
	void closeEvent(QCloseEvent *e);
	void paintEvent(QPaintEvent *paintevent);

private slots:

private:
	void AboutImageWarping();
	void AboutGCL();
	void CreateActions();
	void CreateMenus();
	void CreateToolBars();
	void CreateStatusBar();

private:
	QMenu		*menu_file_;
	QMenu		*menu_edit_;
	QMenu		*menu_help_;
	QMenu		*about_menu_;
	QToolBar	*toolbar_file_;
	QToolBar	*toolbar_edit_;
	QAction		*minidraw_about_action_;
	QAction		*gcl_about_action_;
	QAction		*action_new_;
	QAction		*action_open_;
	QAction		*action_save_;
	QAction		*action_saveas_;
	QAction		*action_restore_;
	QAction		*action_exit_;

	QToolBar	*toolbar_parameterization_;
	QAction		*action_pt_uniform;
	QAction		*action_pt_weightedleastsquares;
	QAction		*action_pt_shapepreserving;
	QAction		*action_pt_asap;
	QAction		*action_pt_arap;
	QAction		*action_pt_slim;
};

#endif // MAINWINDOW_H
