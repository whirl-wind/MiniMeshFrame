#include "multiinputdialog.h"
#include <QBoxLayout>
#include <QPushButton>

MultiInputDialog::MultiInputDialog(int count, QWidget * parent)
	: QDialog(parent), m_GroupCount(count)
{
	QVBoxLayout *layout = new QVBoxLayout;
	QHBoxLayout *subLayout;
	m_Labels = new QLabel*[m_GroupCount];
	m_LineEdits = new QLineEdit*[m_GroupCount];
	//��ƽ���
	for (int i = 0; i < m_GroupCount; i++)
	{
		subLayout = new QHBoxLayout;
		m_LineEdits[i] = new QLineEdit(this);
		m_Labels[i] = new QLabel(this);
		subLayout->addWidget(m_Labels[i]);
		subLayout->addWidget(m_LineEdits[i]);
		layout->addLayout(subLayout);
	}

	m_OKButton = new QPushButton(tr("OK"), this);
	m_CancelButton = new QPushButton(tr("Cancel"), this);
	subLayout = new QHBoxLayout;
	subLayout->addStretch();
	subLayout->addWidget(m_OKButton);
	subLayout->addWidget(m_CancelButton);
	layout->addLayout(subLayout);
	setLayout(layout);

	connect(m_OKButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(m_CancelButton, SIGNAL(clicked()), this, SLOT(accept()));
}

MultiInputDialog::~MultiInputDialog()
{
	delete[] m_Labels;
	delete[] m_LineEdits;
	delete m_OKButton;
	delete m_CancelButton;
	m_OKButton = nullptr;
	m_CancelButton = nullptr;
}

void MultiInputDialog::SetLabelTexts(const QStringList & listText)
{
	for (int i = 0; i < listText.size(); i++)
	{
		if (i >= m_GroupCount)
			break;
		m_Labels[i]->setText(listText.at(i));
	}
}

void MultiInputDialog::SetOneLabelText(int index, const QString & text)
{
	m_Labels[index]->setText(text);
}

QString MultiInputDialog::GetOneText(int index)
{
	return m_LineEdits[index]->text();
}

QStringList MultiInputDialog::GetAllTexts()
{
	QStringList list;
	for (int i = 0; i < m_GroupCount; i++)
	{
		list.push_back(m_LineEdits[i]->text());
	}
	return list;
}

void MultiInputDialog::SetLabelsWidth(int width)
{
	for (int i = 0; i < m_GroupCount; i++)
		m_Labels[i]->setFixedWidth(width);
}

void MultiInputDialog::SetLineEditRegExp(int index, QRegExp regExp)
{
	QValidator *validator = new QRegExpValidator(regExp, this);
	m_LineEdits[index]->setValidator(validator);
}
