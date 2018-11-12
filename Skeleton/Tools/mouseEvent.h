#ifndef MOUSEEVENT_H
#define MOUSEEVENT_H

#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  

using namespace std;
using namespace cv;

class mouseEventTool
{
	

private:
	int rad = 60;
	IplImage* src = 0;
	IplImage* img = 0;
	IplImage* dst = 0;

public:
	mouseEventTool();
	virtual ~mouseEventTool();

	void on_mouse(int event, int x, int y, int flags, void* ustc);

	void set_rad(int r);
	int get_rad();

	void set_src(IplImage* src_in);
	IplImage* get_src();

	void set_img(IplImage* img_in);
	IplImage* get_img();

	void set_dst(IplImage* dst_in);
	IplImage* get_dst();


};
#endif
