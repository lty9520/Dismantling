#ifndef BENCHDELETE_H
#define BENCHDELETE_H

#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <list>
#include <iostream>  
#include <algorithm>  
#include <vector>

using namespace std;
using namespace cv;

class benchDelete
{



public:
	benchDelete();
	virtual ~benchDelete();
	//分支点标记
	void bench_delete(Mat srcimg);

};

#endif
