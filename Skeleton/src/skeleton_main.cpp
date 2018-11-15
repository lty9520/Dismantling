#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <list>
#include <iostream>  
#include <algorithm>  
#include <set>
#include <map>
#include <fstream>
#include <stack>
#include "skeleton_tools.h"
#include "benchDelete.h"
#include "mouseEvent.h"
#include "NeibAreaMarking.h"

using namespace std;
using namespace cv;



int main()
{
	cv::Mat raw = cv::imread("333.jpg", 0);
	cv::Mat binaryImage;
	cv::threshold(raw, binaryImage, 180, 255, CV_THRESH_BINARY_INV);
	cv::imshow("二值化图像", binaryImage * 255);
	cout << "first" << endl;
	skeletonThinning st;
	//thin(binaryImage);
	//st.thinning(binaryImage, binaryImage);
	st.zhang_thinimage_improve(binaryImage);
	Mat skeleton_img = binaryImage.clone() * 255;
	cv::imshow("骨架图像", skeleton_img);
	cout << "second" << endl;

	
	benchDelete bd;
	bd.bench_delete(skeleton_img);

	imshow("bench delete", skeleton_img);
	cout << "third" << endl;

	Mat labelImg;
	int* labelmap[50] = {0};
	nbaMarking nbam;
	//nbam.Two_Pass(skeleton_img, labelImg);
	nbam.Seed_Filling(skeleton_img, labelImg, 8);
	//nbam.bwLabel(skeleton_img, labelImg); 

	cout << "4th" << endl;

	Mat colorLabelImg;
	nbam.LabelColor(labelImg, colorLabelImg);
	imshow("colorImg", colorLabelImg);
	cout << "5th" << endl;

	//mouseEventTool met;
	//IplImage* src = 0;
	//IplImage* img = 0;
	//IplImage* dst = 0;
	//IplImage orgTemp = colorLabelImg;
	//src = cvCloneImage(&orgTemp);
	//img = cvCloneImage(src);
	//dst = cvCreateImage(cvSize(foo * 4, foo * 4), src->depth, src->nChannels);
	//met.set_src(src);
	//met.set_dst(dst);
	//met.set_img(img);
	//cvNamedWindow("img", 1);
	//cvSetMouseCallback("img", on_mouse, 0);
	//cvShowImage("img", img);
	//cvNamedWindow("dst", 1);

	cv::cvtColor(raw, raw, CV_GRAY2BGR);
	for (int row = 0; row < binaryImage.rows; row++)
		for (int col = 0; col < binaryImage.cols; col++)
		{
			if (binaryImage.at<uchar>(row, col) == 255)
			{
				raw.at<cv::Vec3b>(row, col)[0] = 255;	//b
				raw.at<cv::Vec3b>(row, col)[1] = 0;		//r
				raw.at<cv::Vec3b>(row, col)[1] = 0;		//g
			}
		}
	cv::imshow("骨架在原图的位置", raw);
	cv::waitKey(0);
	cvDestroyAllWindows();
	//cvReleaseImage(&src);
	//cvReleaseImage(&img);
	//cvReleaseImage(&dst);
	return 0;
}