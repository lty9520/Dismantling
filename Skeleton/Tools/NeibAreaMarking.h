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
	//返回较小的值
	int smallerone(int num1, int num2);
	


public:
	nbaMarking();
	virtual ~nbaMarking();
	//实现的bwLabel算法
	void bwLabel(const Mat& imgBw, Mat & imgLabeled);
	//种子区域生长方法
	void Seed_Filling(const cv::Mat& binImg, cv::Mat& lableImg, int model_flag);
	//两次扫描方法
	void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg);
	//不同连通域标记颜色
	void LabelColor(const cv::Mat& labelImg, cv::Mat& colorLabelImg);


};


#endif
