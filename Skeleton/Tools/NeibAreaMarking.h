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


};


#endif