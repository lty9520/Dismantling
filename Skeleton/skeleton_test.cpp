#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <list>
#include <iostream>  
#include <algorithm>  
#include <set>
#include <map>
#include <fstream>

using namespace std;
using namespace cv;


IplImage* src = 0;
IplImage* img = 0;
IplImage* dst = 0;
int foo = 60;
void on_mouse(int event, int x, int y, int flags, void* ustc)
{
	if (event == CV_EVENT_MOUSEMOVE && !(flags & CV_EVENT_FLAG_LBUTTON))
	{
		int n = 1;
		cvCopy(src, img);
		CvPoint p0;
		CvPoint p1;
		if (x < foo)
		{
			if (y < foo)
			{
				p0 = cvPoint(0, 0);
				p1 = cvPoint(n * foo, n * foo);
			}
			else if (y > img->height - foo)
			{
				p0 = cvPoint(0, img->height - n * foo);
				p1 = cvPoint(n * foo, img->height);
			}
			else
			{
				p0 = cvPoint(0, y - foo);
				p1 = cvPoint(n * foo, y + foo);
			}
		}
		else if (x > img->width - foo)
		{
			if (y < foo)
			{
				p0 = cvPoint(img->width - n * foo, 0);
				p1 = cvPoint(img->width, n * foo);
			}
			else if (y > img->height - foo)
			{
				p0 = cvPoint(img->width - n * foo, img->height - n * foo);
				p1 = cvPoint(img->width, img->height);
			}
			else
			{
				p0 = cvPoint(img->width - n * foo, y - foo);
				p1 = cvPoint(img->width, y + foo);
			}
		}
		else
		{
			if (y < foo)
			{
				p0 = cvPoint(x - foo, 0);
				p1 = cvPoint(x + foo, n * foo);
			}
			else if (y > img->height - foo)
			{
				p0 = cvPoint(x - foo, img->height - n * foo);
				p1 = cvPoint(x + foo, img->height);
			}
			else
			{
				p0 = cvPoint(x - foo, y - foo);
				p1 = cvPoint(x + foo, y + foo);
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

//获取A0~A5
//************************************
// Method:    GetAi
// FullName:  GetAi
// Access:    public 
// Returns:   set<int>
// Qualifier:
// Parameter: int a[]
// Parameter: int length
//************************************
set<int> GetAi(int a[], int length)
{
	set<int> vec;
	int neighbour[] = { 1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64 };
	for (int i = 0; i < length; i++)
		for (int j = 0; j < 8; j++)
		{
			int sum = 0;
			for (int k = j; k <= j + a[i]; k++)
				sum += neighbour[k];
			vec.insert(sum);
			std::cout << sum << " ";
		}
	std::cout << std::endl;
	return vec;
}

//迭代腐蚀
//************************************
// Method:    erodephase
// FullName:  erodephase
// Access:    public 
// Returns:   bool
// Qualifier:
// Parameter: list<cv::Point> & border
// Parameter: cv::Mat & Input
// Parameter: int neighbour[][3]
// Parameter: const set<int> & A
//************************************
bool erodephase(list<cv::Point> &border, cv::Mat&Input, int neighbour[][3], const set<int>& A)
{
	auto pt = border.begin();
	bool result = false;
	while (pt != border.end())
	{

		int weight = 0;
		for (int j = -1; j <= 1; ++j)
			for (int k = -1; k <= 1; k++)
				weight += neighbour[j + 1][k + 1] * Input.at<uchar>(pt->y + j, pt->x + k);

		if (std::find(A.begin(), A.end(), weight) != A.end())
		{
			Input.at<uchar>(pt->y, pt->x) = 0;
			pt = border.erase(pt);
			result = true;
		}
		else
			++pt;
	}
	return result;
}

//找边界 
//************************************
// Method:    findborder
// FullName:  findborder
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: list<cv::Point2i> & border
// Parameter: const cv::Mat & Input
//************************************
void findborder(list<cv::Point2i>& border, const cv::Mat&Input)
{
	int cnt = 0;
	int rows = Input.rows;
	int cols = Input.cols;
	cv::Mat bordermat = Input.clone();
	for (int row = 1; row < rows - 1; ++row)
		for (int col = 1; col < cols - 1; ++col)
		{
			int weight = 0;
			for (int j = -1; j <= 1; ++j)
				for (int k = -1; k <= 1; k++)
				{
					if (Input.at<uchar>(row + j, col + k) == 1)
						++cnt;
				}
			if (cnt == 9)
				bordermat.at<uchar>(row, col) = 0;
			cnt = 0;
		}

	for (int row = 1; row < rows - 1; ++row)
		for (int col = 1; col < cols - 1; ++col)
		{
			if (bordermat.at<uchar>(row, col) == 1)
				border.push_back(cv::Point2i(col, row));
		}

}

//最后一步，得到骨架
//************************************
// Method:    finalerode
// FullName:  finalerode
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: cv::Mat & Input
// Parameter: int neighbour[][3]
// Parameter: const set<int> & A
//************************************
void finalerode(cv::Mat&Input, int neighbour[][3], const set<int>& A)
{
	int rows = Input.rows;
	int cols = Input.cols;
	for (int m = 1; m < rows - 1; ++m)
		for (int n = 1; n < cols - 1; ++n)
		{
			int weight = 0;
			for (int j = -1; j <= 1; ++j)
				for (int k = -1; k <= 1; k++)
				{
					weight += neighbour[j + 1][k + 1] * Input.at<uchar>(m + j, n + k);
				}

			if (std::find(A.begin(), A.end(), weight) != A.end())
				Input.at<uchar>(m, n) = 0;
		}
}

//************************************
// Method:    thin
// FullName:  thin
// Access:    public 
// Returns:   void
// Qualifier: //Input是二值图像
// Parameter: cv::Mat & Input
//************************************
void thin(cv::Mat &Input) //Input是二值图像
{
	int a0[] = { 1, 2, 3, 4, 5, 6 };
	int a1[] = { 2 };
	int a2[] = { 2, 3 };
	int a3[] = { 2, 3, 4 };
	int a4[] = { 2, 3, 4, 5 };
	int a5[] = { 2, 3, 4, 5, 6 };
	set<int> A0 = GetAi(a0, 6);
	set<int> A1 = GetAi(a1, 1);
	set<int> A2 = GetAi(a2, 2);
	set<int> A3 = GetAi(a3, 3);
	set<int> A4 = GetAi(a4, 4);
	set<int> A5 = GetAi(a5, 5);
	list<cv::Point2i> border;
	bool continue_ = true;
	int neighbour[3][3] = {
		{ 128, 1, 2 },
		{ 64, 0, 4 },
		{ 32, 16, 8 }
	};
	while (continue_)
	{
		//找边界，对原文方法做了小改变，但影响不大。
		continue_ = false;

		findborder(border, Input);//Phase0
		//可以在下面每一步打印结果，看每一步对提取骨架的贡献
		erodephase(border, Input, neighbour, A1);//Phase1
		erodephase(border, Input, neighbour, A2);//Phase2
		erodephase(border, Input, neighbour, A3);//Phase3
		erodephase(border, Input, neighbour, A4);//Phase4
		continue_ = erodephase(border, Input, neighbour, A5);//Phase5
		border.clear();

	}
	finalerode(Input, neighbour, A0);//最后一步

}

/**
* Perform one thinning iteration.
* Normally you wouldn't call this function directly from your code.
*
* Parameters:
* 		im    Binary image with range = [0,1]
* 		iter  0=even, 1=odd
*/
void thinningIteration(cv::Mat& img, int iter)
{
	CV_Assert(img.channels() == 1);
	CV_Assert(img.depth() != sizeof(uchar));
	CV_Assert(img.rows > 3 && img.cols > 3);

	cv::Mat marker = cv::Mat::zeros(img.size(), CV_8UC1);

	int nRows = img.rows;
	int nCols = img.cols;

	if (img.isContinuous()) {
		nCols *= nRows;
		nRows = 1;
	}

	int x, y;
	uchar *pAbove;
	uchar *pCurr;
	uchar *pBelow;
	uchar *nw, *no, *ne;    // north (pAbove)
	uchar *we, *me, *ea;
	uchar *sw, *so, *se;    // south (pBelow)

	uchar *pDst;

	// initialize row pointers
	pAbove = NULL;
	pCurr = img.ptr<uchar>(0);
	pBelow = img.ptr<uchar>(1);

	for (y = 1; y < img.rows - 1; ++y) {
		// shift the rows up by one
		pAbove = pCurr;
		pCurr = pBelow;
		pBelow = img.ptr<uchar>(y + 1);

		pDst = marker.ptr<uchar>(y);

		// initialize col pointers
		no = &(pAbove[0]);
		ne = &(pAbove[1]);
		me = &(pCurr[0]);
		ea = &(pCurr[1]);
		so = &(pBelow[0]);
		se = &(pBelow[1]);

		for (x = 1; x < img.cols - 1; ++x) {
			// shift col pointers left by one (scan left to right)
			nw = no;
			no = ne;
			ne = &(pAbove[x + 1]);
			we = me;
			me = ea;
			ea = &(pCurr[x + 1]);
			sw = so;
			so = se;
			se = &(pBelow[x + 1]);

			int A = (*no == 0 && *ne == 1) + (*ne == 0 && *ea == 1) +
				(*ea == 0 && *se == 1) + (*se == 0 && *so == 1) +
				(*so == 0 && *sw == 1) + (*sw == 0 && *we == 1) +
				(*we == 0 && *nw == 1) + (*nw == 0 && *no == 1);
			int B = *no + *ne + *ea + *se + *so + *sw + *we + *nw;
			int m1 = iter == 0 ? (*no * *ea * *so) : (*no * *ea * *we);
			int m2 = iter == 0 ? (*ea * *so * *we) : (*no * *so * *we);

			if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)//腐蚀四个条件
				pDst[x] = 1;
		}
	}

	img &= ~marker;
}

/**
* Function for thinning the given binary image
*
* Parameters:
* 		src  The source image, binary with range = [0,255]
* 		dst  The destination image
*/
void thinning(const cv::Mat& src, cv::Mat& dst)
{
	dst = src.clone();
	dst /= 255;         // convert to binary image

	cv::Mat prev = cv::Mat::zeros(dst.size(), CV_8UC1);
	cv::Mat diff;

	do {
		thinningIteration(dst, 0);
		thinningIteration(dst, 1);
		cv::absdiff(dst, prev, diff);
		dst.copyTo(prev);
		static int ii = 0;
		std::cout << ++ii << std::endl;
	} while (cv::countNonZero(diff) >0);//迭代终止条件

	dst *= 255;
}

//************************************
// Method:    zhang_thinimage_improve
// FullName:  zhang_thinimage_improve
// Access:    public 
// Returns:   void
// Qualifier: //单通道、二值化后的图像
// Parameter: Mat & srcimage
//************************************
void zhang_thinimage_improve(Mat &srcimage)//单通道、二值化后的图像  
{
	vector<Point> deletelist1;
	int Zhangmude[9];
	int nl = srcimage.rows;
	int nc = srcimage.cols;
	while (true)
	{
		for (int j = 1; j < (nl - 1); j++)
		{
			uchar* data_last = srcimage.ptr<uchar>(j - 1);
			uchar* data = srcimage.ptr<uchar>(j);
			uchar* data_next = srcimage.ptr<uchar>(j + 1);
			for (int i = 1; i < (nc - 1); i++)
			{
				if (data[i] == 255)
				{
					Zhangmude[0] = 1;
					if (data_last[i] == 255) Zhangmude[1] = 1;
					else  Zhangmude[1] = 0;
					if (data_last[i + 1] == 255) Zhangmude[2] = 1;
					else  Zhangmude[2] = 0;
					if (data[i + 1] == 255) Zhangmude[3] = 1;
					else  Zhangmude[3] = 0;
					if (data_next[i + 1] == 255) Zhangmude[4] = 1;
					else  Zhangmude[4] = 0;
					if (data_next[i] == 255) Zhangmude[5] = 1;
					else  Zhangmude[5] = 0;
					if (data_next[i - 1] == 255) Zhangmude[6] = 1;
					else  Zhangmude[6] = 0;
					if (data[i - 1] == 255) Zhangmude[7] = 1;
					else  Zhangmude[7] = 0;
					if (data_last[i - 1] == 255) Zhangmude[8] = 1;
					else  Zhangmude[8] = 0;
					int whitepointtotal = 0;
					for (int k = 1; k < 9; k++)
					{
						//得到1的个数
						whitepointtotal = whitepointtotal + Zhangmude[k];
					}
					if ((whitepointtotal >= 2) && (whitepointtotal <= 6))
					{
						//得到01的个数
						int ap = 0;
						if ((Zhangmude[1] == 0) && (Zhangmude[2] == 1)) ap++;
						if ((Zhangmude[2] == 0) && (Zhangmude[3] == 1)) ap++;
						if ((Zhangmude[3] == 0) && (Zhangmude[4] == 1)) ap++;
						if ((Zhangmude[4] == 0) && (Zhangmude[5] == 1)) ap++;
						if ((Zhangmude[5] == 0) && (Zhangmude[6] == 1)) ap++;
						if ((Zhangmude[6] == 0) && (Zhangmude[7] == 1)) ap++;
						if ((Zhangmude[7] == 0) && (Zhangmude[8] == 1)) ap++;
						if ((Zhangmude[8] == 0) && (Zhangmude[1] == 1)) ap++;
						//计算bp
						int bp = 0;
						bp += Zhangmude[1];
						bp += Zhangmude[2] << 1;
						bp += Zhangmude[3] << 2;
						bp += Zhangmude[4] << 3;
						bp += Zhangmude[5] << 4;
						bp += Zhangmude[6] << 5;
						bp += Zhangmude[7] << 6;
						bp += Zhangmude[8] << 7;

						if (ap == 1 || bp == 65 || bp == 5 || bp == 20 || bp == 80 || bp == 13 || bp == 22 || bp == 52 || bp == 133 || bp == 141 || bp == 54)
						{
							if ((Zhangmude[1] * Zhangmude[7] * Zhangmude[5] == 0) && (Zhangmude[3] * Zhangmude[5] * Zhangmude[7] == 0))
							{
								deletelist1.push_back(Point(i, j));  //满足条件，去除该点
							}
						}
					}
				}
			}
		}
		if (deletelist1.size() == 0) break;
		for (size_t i = 0; i < deletelist1.size(); i++)
		{
			Point tem;
			tem = deletelist1[i];
			uchar* data = srcimage.ptr<uchar>(tem.y);
			data[tem.x] = 0;
		}
		deletelist1.clear();


		for (int j = 1; j < (nl - 1); j++)
		{
			uchar* data_last = srcimage.ptr<uchar>(j - 1);
			uchar* data = srcimage.ptr<uchar>(j);
			uchar* data_next = srcimage.ptr<uchar>(j + 1);
			for (int i = 1; i < (nc - 1); i++)
			{
				if (data[i] == 255)
				{
					Zhangmude[0] = 1;
					if (data_last[i] == 255) Zhangmude[1] = 1;
					else  Zhangmude[1] = 0;
					if (data_last[i + 1] == 255) Zhangmude[2] = 1;
					else  Zhangmude[2] = 0;
					if (data[i + 1] == 255) Zhangmude[3] = 1;
					else  Zhangmude[3] = 0;
					if (data_next[i + 1] == 255) Zhangmude[4] = 1;
					else  Zhangmude[4] = 0;
					if (data_next[i] == 255) Zhangmude[5] = 1;
					else  Zhangmude[5] = 0;
					if (data_next[i - 1] == 255) Zhangmude[6] = 1;
					else  Zhangmude[6] = 0;
					if (data[i - 1] == 255) Zhangmude[7] = 1;
					else  Zhangmude[7] = 0;
					if (data_last[i - 1] == 255) Zhangmude[8] = 1;
					else  Zhangmude[8] = 0;
					int whitepointtotal = 0;
					for (int k = 1; k < 9; k++)
					{
						whitepointtotal = whitepointtotal + Zhangmude[k];
					}
					if ((whitepointtotal >= 2) && (whitepointtotal <= 6))
					{
						int ap = 0;
						if ((Zhangmude[1] == 0) && (Zhangmude[2] == 1)) ap++;
						if ((Zhangmude[2] == 0) && (Zhangmude[3] == 1)) ap++;
						if ((Zhangmude[3] == 0) && (Zhangmude[4] == 1)) ap++;
						if ((Zhangmude[4] == 0) && (Zhangmude[5] == 1)) ap++;
						if ((Zhangmude[5] == 0) && (Zhangmude[6] == 1)) ap++;
						if ((Zhangmude[6] == 0) && (Zhangmude[7] == 1)) ap++;
						if ((Zhangmude[7] == 0) && (Zhangmude[8] == 1)) ap++;
						if ((Zhangmude[8] == 0) && (Zhangmude[1] == 1)) ap++;
						int bp = 0;
						bp += Zhangmude[1];
						bp += Zhangmude[2] << 1;
						bp += Zhangmude[3] << 2;
						bp += Zhangmude[4] << 3;
						bp += Zhangmude[5] << 4;
						bp += Zhangmude[6] << 5;
						bp += Zhangmude[7] << 6;
						bp += Zhangmude[8] << 7;
						if (ap == 1 || bp == 65 || bp == 5 || bp == 20 || bp == 80 || bp == 13 || bp == 22 || bp == 52 || bp == 133 || bp == 141 || bp == 54)
						{
							if ((Zhangmude[1] * Zhangmude[3] * Zhangmude[5] == 0) && (Zhangmude[3] * Zhangmude[1] * Zhangmude[7] == 0))
							{
								deletelist1.push_back(Point(i, j));
							}
						}
					}
				}
			}
		}
		if (deletelist1.size() == 0) break;
		for (size_t i = 0; i < deletelist1.size(); i++)
		{
			Point tem;
			tem = deletelist1[i];
			uchar* data = srcimage.ptr<uchar>(tem.y);
			data[tem.x] = 0;
		}
		deletelist1.clear();
	}
}

void bench_delete(Mat srcimg)
{
	//骨架端点向量
	vector<Point> end_point;
	//骨架分支点向量
	vector<Point> bench_point;
	//既不是分支点也不是端点向量
	vector<Point> other_point;
	//正常骨架连续点向量
	vector<Point> skele_point;
	int neighbourhood[9];
	int row = srcimg.rows;
	int col = srcimg.cols;
	int num = 0;
	for (int j = 1; j < row - 1; j++)
	{
		uchar* data_last = srcimg.ptr<uchar>(j - 1);
		uchar* data = srcimg.ptr<uchar>(j);
		uchar* data_next = srcimg.ptr<uchar>(j + 1);
		for (int i = 1; i < col - 1; i++)
		{
			if (data[i] == 255)
			{
				//						————————————对应关系————————————
				//	p点8邻域图	|	neighbourhood[]		|								data[]	
				//	9	2	3	|		8	1	2		|	data_last[i - 1]		data_last[i]		data_last[i + 1]
				//	8	1	4	|		7	0	3		|	  data[i - 1]				data[i]				data[i + 1]
				//	7	6	5	|		6	5	4		|	data_next[i - 1]		data_next[i]		data_next[i + 1]
				//
				neighbourhood[0] = 1;
				
				if (data_last[i] == 255)
					neighbourhood[1] = 1;
				else
					neighbourhood[1] = 0;

				if (data_last[i + 1] == 255)
					neighbourhood[2] = 1;
				else
					neighbourhood[2] = 0;

				if (data[i + 1] == 255)
					neighbourhood[3] = 1;
				else
					neighbourhood[3] = 0;

				if (data_next[i + 1] == 255)
					neighbourhood[4] = 1;
				else
					neighbourhood[4] = 0;

				if (data_next[i] == 255)
					neighbourhood[5] = 1;
				else
					neighbourhood[5] = 0;

				if (data_next[i - 1] == 255)
					neighbourhood[6] = 1;
				else
					neighbourhood[6] = 0;

				if (data[i - 1] == 255)
					neighbourhood[7] = 1;
				else
					neighbourhood[7] = 0;

				if (data_last [i - 1] == 255)
					neighbourhood[8] = 1;
				else
					neighbourhood[8] = 0;

				int neighbour_sum = 0;
				for (int k = 1; k < 9; k++)
				{
					neighbour_sum += neighbourhood[k];
				}

				if (neighbourhood[0] == 1 && neighbour_sum < 2)
					end_point.push_back(Point(i, j));
				else if (neighbourhood[0] == 1 && neighbour_sum > 2)
				{
					bench_point.push_back(Point(i, j));
					data[i] = 1;
					num++;
				}
				else if (neighbourhood[0] == 1 && neighbour_sum == 2 &&
					(neighbourhood[1] + neighbourhood[2] == 2 ||
					neighbourhood[2] + neighbourhood[3] == 2 ||
					neighbourhood[3] + neighbourhood[4] == 2 ||
					neighbourhood[4] + neighbourhood[5] == 2 ||
					neighbourhood[5] + neighbourhood[6] == 2 ||
					neighbourhood[6] + neighbourhood[7] == 2 ||
					neighbourhood[7] + neighbourhood[8] == 2
					))
					other_point.push_back(Point(i, j));
				else
					skele_point.push_back(Point(i, j));
			}
		}
	}
	cout << "分支点个数:" << num << endl;
	
}


void bwLabel(const Mat& imgBw, Mat & imgLabeled)
{
	// 对图像周围扩充一格
	Mat imgClone = Mat(imgBw.rows + 1, imgBw.cols + 1, imgBw.type(), Scalar(0));
	imgBw.copyTo(imgClone(Rect(1, 1, imgBw.cols, imgBw.rows)));

	imgLabeled.create(imgClone.size(), imgClone.type());
	imgLabeled.setTo(Scalar::all(0));

	vector<vector<Point>> contours;
	vector<Vec4i> hierarchy;
	findContours(imgClone, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);

	vector<int> contoursLabel(contours.size(), 0);
	int numlab = 1;
	// 标记外围轮廓
	for (vector<vector<Point>>::size_type i = 0; i < contours.size(); i++)
	{
		if (hierarchy[i][3] >= 0) // 有父轮廓
		{
			continue;
		}
		for (vector<Point>::size_type k = 0; k != contours[i].size(); k++)
		{
			imgLabeled.at<uchar>(contours[i][k].y, contours[i][k].x) = numlab;
		}
		contoursLabel[i] = numlab++;
	}
	// 标记内轮廓
	for (vector<vector<Point>>::size_type i = 0; i < contours.size(); i++)
	{
		if (hierarchy[i][3] < 0)
		{
			continue;
		}
		for (vector<Point>::size_type k = 0; k != contours[i].size(); k++)
		{
			imgLabeled.at<uchar>(contours[i][k].y, contours[i][k].x) = contoursLabel[hierarchy[i][3]];
		}
	}
	// 非轮廓像素的标记
	for (int i = 0; i < imgLabeled.rows; i++)
	{
		for (int j = 0; j < imgLabeled.cols; j++)
		{
			if (imgClone.at<uchar>(i, j) != 0 && imgLabeled.at<uchar>(i, j) == 0)
			{
				imgLabeled.at<uchar>(i, j) = imgLabeled.at<uchar>(i, j - 1);
			}
		}
	}
	imgLabeled = imgLabeled(Rect(1, 1, imgBw.cols, imgBw.rows)).clone(); // 将边界裁剪掉1像素
}


int smallerone(int num1, int num2)
{
	if (num1 < num2)
		return num1;
	else
		return num2;
}


void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg)    //两遍扫描法
{
	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	//IplImage binImgTemp = binImg;
	//src = cvCloneImage(&binImgTemp);
	//IplImage *pGrayImage_8U = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 1);
	//IplImage *pGrayImage_32S = cvCreateImage(cvGetSize(src), IPL_DEPTH_32S, 1);
	//
	//cvConvertScale(pGrayImage_8U, pGrayImage_32S);  //8U转32S
	// 第一个通路

	lableImg.release();
	binImg.convertTo(lableImg, CV_8UC1);

	//lableImg = cvarrToMat(pGrayImage_32S);

	int label = 1;
	std::vector<int> labelSet;
	labelSet.push_back(0);
	labelSet.push_back(1);

	//uchar *temp = lableImg.ptr<uchar>(26);//获取第五行的首地址
	//int a = temp[248];//获取第五行第五列的像素值并幅值给a。
	//cout << a << endl;//输出像素值
	//ofstream fout("1.txt");
	//
	//for (int nrow = 0; nrow < binImg.rows; nrow++)
	//{
	//	for (int ncol = 0; ncol < binImg.cols; ncol++)
	//	{
	//		uchar val = binImg.at<uchar>(nrow, ncol);
	//		if (val == 255)
	//		{
	//			val = 1;
	//		}
	//		else
	//		{
	//			val = 0;
	//		}
	//
	//		fout << int(val);		//将像素值转化为整数后写入文件
	//		//cout << int(val);
	//	}
	//	//cout << endl;			//打印一行后进行换行
	//	fout << endl;
	//}
	//
	//fout.close();

	//ofstream fout("3.txt");
	//for (int i = 0; i < lableImg.rows; i++)
	//{
	//	uchar* temp = lableImg.ptr<uchar>(i);
	//	for (int j = 0; j < lableImg.cols; j++)
	//	{
	//		int a = temp[j];
	//		fout << a;
	//		fout << "\t";
	//	}
	//	fout << "\n";
	//}
	//fout.close();
	

	int rows = binImg.rows - 1;
	int cols = binImg.cols - 1;
	for (int i = 1; i < rows; i++)
	{
		int* data_preRow = lableImg.ptr<int>(i - 1);
		int* data_curRow = lableImg.ptr<int>(i);
		for (int j = 1; j < cols; j++)
		{
			if (data_curRow[j] == 1)
			{
				std::vector<int> neighborLabels;
				neighborLabels.reserve(2);
				int leftPixel = data_curRow[j - 1];
				int upPixel = data_preRow[j];
				int leftupPixel = data_preRow[j - 1];
				int rightupPixel = data_preRow[j + 1];
				
				if (leftPixel > 1 && rightupPixel > 1)
				{
					if (leftPixel < rightupPixel)
					{
						neighborLabels.push_back(leftPixel);
						rightupPixel = leftPixel;
					} 
					else
					{
						neighborLabels.push_back(rightupPixel);
						leftPixel = rightupPixel;
					}
				}
				if (leftupPixel > 1 && rightupPixel > 1)
				{
					if (leftupPixel < rightupPixel)
					{
						neighborLabels.push_back(leftupPixel);
						rightupPixel = leftupPixel;
					}
					else
					{
						neighborLabels.push_back(rightupPixel);
						leftupPixel = rightupPixel;
					}
				}

				if (neighborLabels.empty())
				{
					labelSet.push_back(++label);  // 不连通，标签+1
					data_curRow[j] = label;
					labelSet[label] = label;
				}
				else
				{
					if (leftPixel > 1)
					{
						neighborLabels.push_back(leftPixel);
					} 
					else if (leftupPixel > 1)
					{
						neighborLabels.push_back(leftupPixel);
					}
					else if (upPixel > 1)
					{
						neighborLabels.push_back(upPixel);
					}
					else if (rightupPixel > 1)
					{
						neighborLabels.push_back(rightupPixel);
					}

					std::sort(neighborLabels.begin(), neighborLabels.end());
					int smallestLabel = neighborLabels[0];


					// 保存最小等价表
					for (size_t k = 1; k < neighborLabels.size(); k++)
					{
						int tempLabel = neighborLabels[k];
						int& oldSmallestLabel = labelSet[tempLabel];
						if (oldSmallestLabel > smallestLabel)
						{
							labelSet[oldSmallestLabel] = smallestLabel;
							oldSmallestLabel = smallestLabel;
						}
						else if (oldSmallestLabel < smallestLabel)
						{
							labelSet[smallestLabel] = oldSmallestLabel;
						}
					}
				}
			}
		}
	}

	// 更新等价对列表
	// 将最小标号给重复区域
	for (size_t i = 2; i < labelSet.size(); i++)
	{
		int curLabel = labelSet[i];
		int preLabel = labelSet[curLabel];
		while (preLabel != curLabel)
		{
			curLabel = preLabel;
			preLabel = labelSet[preLabel];
		}
		labelSet[i] = curLabel;
	};

	//uchar *temp1 = lableImg.ptr<uchar>(78);//获取第五行的首地址
	//int a1 = temp1[204];//获取第五行第五列的像素值并幅值给a。
	//cout << a1 << endl;//输出像素值

	ofstream fout("3.txt");
	for (int i = 0; i < lableImg.rows; i++)
	{
		uchar* temp = lableImg.ptr<uchar>(i);
		for (int j = 0; j < lableImg.cols; j++)
		{
			int a = temp[j];
			fout << a;
			fout << "\t";
		}
		fout << "\n";
	}
	fout.close();


	for (int i = 0; i < rows; i++)
	{
		int* data = lableImg.ptr<int>(i);
		for (int j = 0; j < cols; j++)
		{
			int pixelLabel = data[j];
			pixelLabel = labelSet[pixelLabel];
		}
	}
}


/*

void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg)    //两遍扫描法
{
	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	// 第一个通路

	lableImg.release();
	binImg.convertTo(lableImg, CV_32SC1);

	int label = 1;
	std::vector<int> labelSet;
	labelSet.push_back(0);
	labelSet.push_back(1);

	int rows = binImg.rows - 1;
	int cols = binImg.cols - 1;
	for (int i = 1; i < rows; i++)
	{
		int* data_preRow = lableImg.ptr<int>(i - 1);
		int* data_curRow = lableImg.ptr<int>(i);
		for (int j = 1; j < cols; j++)
		{
			if (data_curRow[j] == 1)
			{
				std::vector<int> neighborLabels;
				neighborLabels.reserve(2);
				int leftPixel = data_curRow[j - 1];
				int upPixel = data_preRow[j];
				if (leftPixel > 1)
				{
					neighborLabels.push_back(leftPixel);
				}
				if (upPixel > 1)
				{
					neighborLabels.push_back(upPixel);
				}

				if (neighborLabels.empty())
				{
					labelSet.push_back(++label);  // 不连通，标签+1
					data_curRow[j] = label;
					labelSet[label] = label;
				}
				else
				{
					std::sort(neighborLabels.begin(), neighborLabels.end());
					int smallestLabel = neighborLabels[0];
					data_curRow[j] = smallestLabel;

					// 保存最小等价表
					for (size_t k = 1; k < neighborLabels.size(); k++)
					{
						int tempLabel = neighborLabels[k];
						int& oldSmallestLabel = labelSet[tempLabel];
						if (oldSmallestLabel > smallestLabel)
						{
							labelSet[oldSmallestLabel] = smallestLabel;
							oldSmallestLabel = smallestLabel;
						}
						else if (oldSmallestLabel < smallestLabel)
						{
							labelSet[smallestLabel] = oldSmallestLabel;
						}
					}
				}
			}
		}
	}

	// 更新等价对列表
	// 将最小标号给重复区域
	for (size_t i = 2; i < labelSet.size(); i++)
	{
		int curLabel = labelSet[i];
		int preLabel = labelSet[curLabel];
		while (preLabel != curLabel)
		{
			curLabel = preLabel;
			preLabel = labelSet[preLabel];
		}
		labelSet[i] = curLabel;
	};

	for (int i = 0; i < rows; i++)
	{
		int* data = lableImg.ptr<int>(i);
		for (int j = 0; j < cols; j++)
		{
			int& pixelLabel = data[j];
			pixelLabel = labelSet[pixelLabel];
		}
	}
}

*/

//彩色显示
cv::Scalar GetRandomColor()
{
	uchar r = 255 * (rand() / (1.0 + RAND_MAX));
	uchar g = 255 * (rand() / (1.0 + RAND_MAX));
	uchar b = 255 * (rand() / (1.0 + RAND_MAX));
	return cv::Scalar(b, g, r);
}


void LabelColor(const cv::Mat& labelImg, cv::Mat& colorLabelImg)
{
	if (labelImg.empty() ||
		labelImg.type() != CV_8UC1)
	{
		return;
	}

	std::map<int, cv::Scalar> colors;

	int rows = labelImg.rows;
	int cols = labelImg.cols;

	colorLabelImg.release();
	colorLabelImg.create(rows, cols, CV_8UC3);
	colorLabelImg = cv::Scalar::all(0);

	for (int i = 0; i < rows; i++)
	{
		const int* data_src = (int*)labelImg.ptr<int>(i);
		uchar* data_dst = colorLabelImg.ptr<uchar>(i);
		for (int j = 0; j < cols; j++)
		{
			int pixelValue = data_src[j];
			if (pixelValue > 1)
			{
				if (colors.count(pixelValue) <= 0)
				{
					colors[pixelValue] = GetRandomColor();
				}

				cv::Scalar color = colors[pixelValue];
				*data_dst++ = color[0];
				*data_dst++ = color[1];
				*data_dst++ = color[2];
			}
			else
			{
				data_dst++;
				data_dst++;
				data_dst++;
			}
		}
	}
}


int main()
{
	cv::Mat raw = cv::imread("333.jpg", 0);
	cv::Mat binaryImage;
	cv::threshold(raw, binaryImage, 180, 255, CV_THRESH_BINARY_INV);
	cv::imshow("二值化图像", binaryImage * 255);
	cout << "first" << endl;
	//thin(binaryImage);
	//thinning(binaryImage, binaryImage);
	zhang_thinimage_improve(binaryImage);
	Mat skeleton_img = binaryImage.clone() * 255;
	cv::imshow("骨架图像", skeleton_img);
	cout << "second" << endl;

	

	bench_delete(skeleton_img);

	imshow("bench delete", skeleton_img);
	cout << "third" << endl;

	Mat labelImg;
	//Two_Pass(skeleton_img, labelImg);
	bwLabel(skeleton_img, labelImg); 

	cout << "4th" << endl;

	Mat colorLabelImg;
	LabelColor(labelImg, colorLabelImg);
	imshow("lableImg", labelImg);
	imshow("colorImg", colorLabelImg);
	cout << "5th" << endl;

	//IplImage orgTemp = colorLabelImg;
	//src = cvCloneImage(&orgTemp);
	//img = cvCloneImage(src);
	//dst = cvCreateImage(cvSize(foo * 4, foo * 4), src->depth, src->nChannels);
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
	cvReleaseImage(&src);
	cvReleaseImage(&img);
	cvReleaseImage(&dst);
	return 0;
}