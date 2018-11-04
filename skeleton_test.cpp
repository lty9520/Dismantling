#include <opencv2/imgproc/imgproc.hpp>  
#include <opencv2/objdetect/objdetect.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include<list>
#include<iostream>  
#include<algorithm>  
#include<set>
using namespace std;

set<int> GetAi(int a[], int length)//获取A0~A5
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
int main()
{
	cv::Mat raw = cv::imread("222.png", 0);
	cv::Mat binaryImage;
	cv::threshold(raw, binaryImage, 180, 1, CV_THRESH_BINARY_INV);
	cv::imshow("二值化图像", binaryImage * 255);
	thin(binaryImage);
	//thinning(binaryImage, binaryImage);
	cv::imshow("骨架图像", binaryImage * 255);

	cv::cvtColor(raw, raw, CV_GRAY2BGR);
	for (int row = 0; row < binaryImage.rows; row++)
		for (int col = 0; col < binaryImage.cols; col++)
		{
			if (binaryImage.at<uchar>(row, col) == 1)
			{
				raw.at<cv::Vec3b>(row, col)[0] = 0;
				raw.at<cv::Vec3b>(row, col)[1] = 255;
				raw.at<cv::Vec3b>(row, col)[1] = 0;
			}
		}
	cv::imshow("骨架在原图的位置", raw);
	cv::waitKey(0);
	return 0;
}