#include "..\Tools\skeleton_tools.h"


skeletonThinning::skeletonThinning()
{

}
skeletonThinning::~skeletonThinning()
{

}
//获取A0~A5
set<int> skeletonThinning::GetAi(int a[], int length)
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
bool skeletonThinning::erodephase(list<cv::Point> &border, cv::Mat&Input, int neighbour[][3], const set<int>& A)
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
void skeletonThinning::findborder(list<cv::Point2i>& border, const cv::Mat&Input)
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
void skeletonThinning::finalerode(cv::Mat&Input, int neighbour[][3], const set<int>& A)
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
//细化操作
void skeletonThinning::thin(cv::Mat &Input) //Input是二值图像
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
//Perform one thinning iteration.
void skeletonThinning::thinningIteration(cv::Mat& img, int iter)
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
//Function for thinning the given binary image
void skeletonThinning::thinning(const cv::Mat& src, cv::Mat& dst)
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
	} while (cv::countNonZero(diff) > 0);//迭代终止条件

	dst *= 255;
}
////单通道、二值化后的图像细化改进方法
void skeletonThinning::zhang_thinimage_improve(Mat &srcimage)//单通道、二值化后的图像  
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
