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
	//分支长度阈值（过长分支剪掉阈值）
	int lenThres_long = 20;
	//分支长度阈值（过短分支剪掉阈值）
	int lenThres_short = 20;
	//过长阈值分支删除
	void bench_cut_long(vector<vector<int> > labels, Mat labelImg, int lenThres_long, int rows, int cols);
	//过短阈值分支删除
	void bench_cut_short(vector<vector<int> > labels, Mat labelImg, int lenThres_long, int rows, int cols);
	


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
	//设置筛选过长分支长度阈值
	void setLenThres_long(int lenThres_long);
	//得到筛选过长分支长度阈值
	int getLenThres_long();
	//设置筛选过短分支长度阈值
	void setLenThres_short(int lenThres_short);
	//得到筛选过短分支长度阈值
	int getLenThres_short();


};


#endif
