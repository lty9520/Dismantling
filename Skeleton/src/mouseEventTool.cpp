#include "..\Tools\mouseEvent.h"

mouseEventTool::mouseEventTool()
{

}

mouseEventTool::~mouseEventTool()
{

}

void mouseEventTool::on_mouse(int event, int x, int y, int flags, void* ustc)
{
	if (event == CV_EVENT_MOUSEMOVE && !(flags & CV_EVENT_FLAG_LBUTTON))
	{
		int n = 1;
		cvCopy(src, img);
		CvPoint p0;
		CvPoint p1;
		if (x < rad)
		{
			if (y < rad)
			{
				p0 = cvPoint(0, 0);
				p1 = cvPoint(n * rad, n * rad);
			}
			else if (y > img->height - rad)
			{
				p0 = cvPoint(0, img->height - n * rad);
				p1 = cvPoint(n * rad, img->height);
			}
			else
			{
				p0 = cvPoint(0, y - rad);
				p1 = cvPoint(n * rad, y + rad);
			}
		}
		else if (x > img->width - rad)
		{
			if (y < rad)
			{
				p0 = cvPoint(img->width - n * rad, 0);
				p1 = cvPoint(img->width, n * rad);
			}
			else if (y > img->height - rad)
			{
				p0 = cvPoint(img->width - n * rad, img->height - n * rad);
				p1 = cvPoint(img->width, img->height);
			}
			else
			{
				p0 = cvPoint(img->width - n * rad, y - rad);
				p1 = cvPoint(img->width, y + rad);
			}
		}
		else
		{
			if (y < rad)
			{
				p0 = cvPoint(x - rad, 0);
				p1 = cvPoint(x + rad, n * rad);
			}
			else if (y > img->height - rad)
			{
				p0 = cvPoint(x - rad, img->height - n * rad);
				p1 = cvPoint(x + rad, img->height);
			}
			else
			{
				p0 = cvPoint(x - rad, y - rad);
				p1 = cvPoint(x + rad, y + rad);
			}
		}
		cvRectangle(img, p0, p1, CV_RGB(0, 255, 0));
		cvSetImageROI(src, cvRect(p0.x, p0.y, p1.x - p0.x, p1.y - p0.y));
		cvResize(src, dst);
		cvResetImageROI(src);
		cvShowImage("img", img);
		cvShowImage("dst", dst);
	}
}

void mouseEventTool::set_rad(int r)
{
	rad = r;
}

int mouseEventTool::get_rad()
{
	return rad;
}
void mouseEventTool::set_src(IplImage* src_in)
{
	src = src_in;
}

IplImage* mouseEventTool::get_src()
{
	return src;
}

void mouseEventTool::set_img(IplImage* img_in)
{
	img = img_in;
}

IplImage* mouseEventTool::get_img()
{
	return img;
}

void mouseEventTool::set_dst(IplImage* dst_in)
{
	dst = dst_in;
}

IplImage* mouseEventTool::get_dst()
{
	return dst;
}