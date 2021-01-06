#include <qevent.h>
#include <qpixmap.h>
#include <qimage.h>
#include <qpainter.h>
#include "colorbar.h"

ColorBar::ColorBar(Qt::Orientation o, QWidget *parent) :
	QWidget(parent),
	d_orientation(o)
{
/*
#ifndef QT_NO_CURSOR
	setCursor(Qt::PointingHandCursor);
#endif
*/
	setbackgroundtrans();
}

void ColorBar::setOrientation(Qt::Orientation o)
{
	d_orientation = o;
	update();
}


void ColorBar::setRange(double z_min, double z_max)
{
	if (z_min > z_max) {
		min = z_max;
		max = z_min;
	}
	else {
		min = z_min;
		max = z_max;
	}
}

QColor ColorBar::getColor(double x)
{
	double colorBarLength = 343.0;//设置颜色条的长度
	QColor color;
	double tempLength = colorBarLength / 4;

	double i = (x - min) / (max - min) * colorBarLength;

	if (i < 0) {
		color.setRgbF(0, 0, 0);
	}
	else if (i < tempLength / 2) {
		color.setRgbF(0, 0, (tempLength / 2 + i) / tempLength);
	}
	else if (i < tempLength / 2 + tempLength) {
		color.setRgbF(0, (i - tempLength / 2) / tempLength, 1);
	}
	else if (i < tempLength / 2 + 2 * tempLength) {
		color.setRgbF((i - tempLength - tempLength / 2) / tempLength, 1, (tempLength * 2 + tempLength / 2 - i) / tempLength);
	}
	else if (i < tempLength / 2 + 3 * tempLength) {
		color.setRgbF(1, (tempLength * 3 + tempLength / 2 - i) / tempLength, 0);
	}
	else if (i <= colorBarLength) {
		color.setRgbF((colorBarLength - i + tempLength / 2) / (tempLength), 0, 0);
	}
	else {
		color.setRgbF(1.0, 1.0, 1.0);
	}

	return color;
}

void ColorBar::setbackground(const QColor & color)
{
	c_background = color;
}

void ColorBar::setbackgroundtrans()
{
	c_background.setRgbF(0, 0, 0, 0);
}

void ColorBar::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	QColor color;
	QRect section;
	float colorBarLength = 343.0;//设置颜色条的长度
	int text_width = 30;

	section.setRect(0, 0, text_width + 20 + 5, colorBarLength + 50 + 5);
	painter.fillRect(section, c_background);
	/*
	//------设置为gray颜色条---------//
	for (int i = 0; i <= colorBarLength; i++)// gray
	{
		//color.setRgbF(i/colorBarLength,i/colorBarLength,i/colorBarLength);//也可以使用这种方法
		color.setHsv(0, 0, (colorBarLength - i) / colorBarLength * 255);
		section.setRect(150, 50 + i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	*/
	//------设置为jet颜色条---------//
	float tempLength = colorBarLength / 4;
	for (int i = 0; i < tempLength / 2; i++)// jet
	{
		color.setRgbF(0, 0, (tempLength / 2 + i) / tempLength);
		section.setRect(text_width, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	for (int i = tempLength / 2 + 1; i < tempLength / 2 + tempLength; i++)// jet
	{
		color.setRgbF(0, (i - tempLength / 2) / tempLength, 1);
		section.setRect(text_width, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	for (int i = tempLength / 2 + tempLength + 1; i < tempLength / 2 + 2 * tempLength; i++)// jet
	{
		color.setRgbF((i - tempLength - tempLength / 2) / tempLength, 1, (tempLength * 2 + tempLength / 2 - i) / tempLength);
		section.setRect(text_width, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	for (int i = tempLength / 2 + 2 * tempLength + 1; i < tempLength / 2 + 3 * tempLength; i++)// jet
	{
		color.setRgbF(1, (tempLength * 3 + tempLength / 2 - i) / tempLength, 0);
		section.setRect(text_width, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	for (int i = tempLength / 2 + 3 * tempLength + 1; i < colorBarLength; i++)// jet
	{
		color.setRgbF((colorBarLength - i + tempLength / 2) / (tempLength), 0, 0);
		section.setRect(text_width, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	/*
	//------设置为hsv颜色条---------//
	for (int i = 0; i <= colorBarLength; i++)// hsv
	{
		color.setHsvF(i / colorBarLength, 1, 1);
		section.setRect(250, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	*/
	/*
	//------设置为hot颜色条---------//
	tempLength = colorBarLength / 2.5;
	for (int i = 0; i < tempLength / 2; i++)// hot
	{
		color.setRgbF((tempLength / 2 + i) / tempLength, 0, 0);
		section.setRect(300, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	for (int i = tempLength / 2 + 1; i < tempLength / 2 + tempLength; i++)// hot
	{
		color.setRgbF(1, (i - tempLength / 2) / tempLength, 0);
		section.setRect(300, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}

	for (int i = tempLength / 2 + tempLength + 1; i < colorBarLength; i++)// hot
	{
		color.setRgbF(1, 1, (i - tempLength / 2 - tempLength) / (colorBarLength - tempLength / 2 - tempLength + 20));
		section.setRect(300, colorBarLength + 50 - i * 1, 20, 1);
		painter.fillRect(section, color);
	}
	*/
	//---------设置边框--------------//

	painter.drawRect(text_width, 50, 20, colorBarLength);
	painter.setFont(QFont(QString::fromLocal8Bit("宋体"), 10, -1, false));
	painter.drawText(text_width - 5, 40, QStringLiteral("Jet"));

	std::string str;
	QString qstr;

	str = std::to_string(float(max));
	qstr = QString::fromStdString(str);
	painter.drawText(0,60, qstr);

	str = std::to_string(float(min));
	qstr = QString::fromStdString(str);
	painter.drawText(0, colorBarLength + 50, qstr);
}
