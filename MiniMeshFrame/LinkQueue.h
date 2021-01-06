#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>
#include<math.h>

#ifndef _LINKQUEUE_H_
#define _LINKQUEUE_H_

#define OK 1
#define ERROR -1
#define TRUE 1
#define FALSE 0
#define Abort(msg)   { printf("aborting:%s\n",msg);exit(0);}

typedef int Status;

template<class QElemType>
class LinkQueue{
private:
	typedef struct QNode {
		QElemType data;
		struct QNode *next;
	}QNode, *QueuePtr;

	QueuePtr front;
	QueuePtr rear;
public:
	LinkQueue();
	~LinkQueue();
	Status InitQueue();
	Status QueueEmpty();
	Status ClearQueue();
	Status DeQueue(QElemType &e);
	Status EnQueue(QElemType e);
};

template<class QElemType>
inline LinkQueue<QElemType>::LinkQueue()
{
	InitQueue();
}

template<class QElemType>
inline LinkQueue<QElemType>::~LinkQueue()
{
}

template<class QElemType>
inline Status LinkQueue<QElemType>::InitQueue()
{
	front = (QNode *)malloc(sizeof(QNode));
	if (!front) Abort("malloc failed.");
	rear = front;
	return OK;
}

template<class QElemType>
inline Status LinkQueue<QElemType>::QueueEmpty(){
	if(front==rear)return TRUE;
	else return FALSE;
} 

template<class QElemType>
inline Status LinkQueue<QElemType>::ClearQueue(){
	if(QueueEmpty()) return OK;
	QueuePtr p;
	p=front;
	front=front->next; 
	free(p);
	ClearQueue();
}

template<class QElemType>
inline Status LinkQueue<QElemType>::DeQueue(QElemType &e){
	if(QueueEmpty()) return ERROR;
	QueuePtr p;
	p=front;
	front=front->next; 
	e=p->data; 
	free(p);
	return OK;
}

template<class QElemType>
inline Status LinkQueue<QElemType>::EnQueue(QElemType e){
	rear->data=e;
	rear->next=(QNode *)malloc(sizeof(QNode));
	rear=rear->next;
    if(!rear) Abort("malloc failed.");
    return OK;
}

#endif

