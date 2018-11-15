#ifndef NEIBAREAMARKING_H
#define NEIBAREAMARKING_H

#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  

#include <list>
#include <iostream>  
#include <algorithm>  
#include <map>
#include <stack>

using namespace std;
using namespace cv;

class nbaMarking
{
private:
	//���ؽ�С��ֵ
	int smallerone(int num1, int num2);
	//��֧������ֵ��������֧������ֵ��
	int lenThres_long = 20;
	//��֧������ֵ�����̷�֧������ֵ��
	int lenThres_short = 20;
	//������ֵ��֧ɾ��
	void bench_cut_long(vector<vector<int> > labels, Mat labelImg, int lenThres_long, int rows, int cols);
	//������ֵ��֧ɾ��
	void bench_cut_short(vector<vector<int> > labels, Mat labelImg, int lenThres_long, int rows, int cols);
	


public:
	nbaMarking();
	virtual ~nbaMarking();
	//ʵ�ֵ�bwLabel�㷨
	void bwLabel(const Mat& imgBw, Mat & imgLabeled);
	//����������������
	void Seed_Filling(const cv::Mat& binImg, cv::Mat& lableImg, int model_flag);
	//����ɨ�跽��
	void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg);
	//��ͬ��ͨ������ɫ
	void LabelColor(const cv::Mat& labelImg, cv::Mat& colorLabelImg);
	//����ɸѡ������֧������ֵ
	void setLenThres_long(int lenThres_long);
	//�õ�ɸѡ������֧������ֵ
	int getLenThres_long();
	//����ɸѡ���̷�֧������ֵ
	void setLenThres_short(int lenThres_short);
	//�õ�ɸѡ���̷�֧������ֵ
	int getLenThres_short();


};


#endif
