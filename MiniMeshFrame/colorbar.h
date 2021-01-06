#ifndef COLORBAR_H_
#define COLORBAR_H_

#include <qwidget.h>

class ColorBar : public QWidget
{
	Q_OBJECT

public:
	ColorBar(Qt::Orientation = Qt::Horizontal, QWidget * = NULL);

	virtual void setOrientation(Qt::Orientation);
	Qt::Orientation orientation() const { return d_orientation; }

	void setRange(double z_min, double z_max);
	QColor getColor(double x);
	void setbackground(const QColor &color);
	void setbackgroundtrans();

Q_SIGNALS:
	void selected(const QColor &);

protected:
	virtual void paintEvent(QPaintEvent *);

private:
	Qt::Orientation d_orientation;
	double min;
	double max;
	QColor c_background;
};

#endif // !COLORBAR_H_