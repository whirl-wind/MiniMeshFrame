#ifndef MULTITNPUTDIALOG_H
#define MULTITNPUTDIALOG_H

#include <QtGui>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>

class MultiInputDialog : public QDialog
{
	Q_OBJECT

private:
	const int m_GroupCount;
	QLabel **m_Labels;
	QLineEdit **m_LineEdits;
	QPushButton *m_OKButton;
	QPushButton *m_CancelButton;

public:
	MultiInputDialog(int count, QWidget *parent = 0);
	virtual ~MultiInputDialog();
	void SetLabelTexts(const QStringList &listText);
	void SetOneLabelText(int index, const QString &text);
	QString GetOneText(int index);
	QStringList GetAllTexts();
	//Ϊ���ö����������������Щ�������ҹ̶���QLabel�Ŀ��
	void SetLabelsWidth(int width);
	//ʹ��������ʽ������������ַ�
	void SetLineEditRegExp(int index, QRegExp regExp);

	//�������Ҫ��д����������
	virtual void accept() { QDialog::accept(); }
	virtual void reject() { QDialog::reject(); }

};

#endif // !MULTITNPUTDIALOG_H