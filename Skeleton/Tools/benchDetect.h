#ifndef BENCHDETECT_H
#define BENCHDETECT_H

#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <list>
#include <iostream>  
#include <algorithm>  
#include <vector>

using namespace std;
using namespace cv;

class benchDetect
{



public:
	benchDetect();
	virtual ~benchDetect();
	//��֧����
	void bench_detect(Mat srcimg);

};

#endif
